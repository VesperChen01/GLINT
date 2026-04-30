#!/usr/bin/env python3
"""
Word to Markdown Converter using python-docx
Alternative to pandoc for converting Word documents to Markdown
"""

import os
import sys
import re
from pathlib import Path

def check_dependencies():
    """Check if required dependencies are installed"""
    try:
        import docx
        print("✅ python-docx found!")
        return True
    except ImportError:
        print("❌ python-docx not found.")
        print("Install with: pip install python-docx")
        return False

def convert_word_to_markdown(input_path, output_path):
    """
    Convert Word document to Markdown format

    Args:
        input_path: Path to input Word document
        output_path: Path to output Markdown file
    """
    try:
        from docx import Document

        print(f"Reading: {input_path}")
        doc = Document(input_path)

        markdown_lines = []
        markdown_lines.append("# GLINT Tutorial\n")
        markdown_lines.append("*Converted from Word document*\n\n")

        for para in doc.paragraphs:
            text = para.text.strip()
            if not text:
                markdown_lines.append("\n")
                continue

            # Handle heading levels
            if para.style.name.startswith('Heading 1'):
                markdown_lines.append(f"# {text}\n\n")
            elif para.style.name.startswith('Heading 2'):
                markdown_lines.append(f"## {text}\n\n")
            elif para.style.name.startswith('Heading 3'):
                markdown_lines.append(f"### {text}\n\n")
            elif para.style.name.startswith('Heading 4'):
                markdown_lines.append(f"#### {text}\n\n")
            else:
                # Regular paragraph - preserve basic formatting
                formatted_text = text

                # Handle bold text (simple heuristic)
                # Note: python-docx doesn't preserve formatting well, so this is basic
                markdown_lines.append(f"{formatted_text}\n\n")

        # Handle tables
        for table in doc.tables:
            markdown_lines.append("\n| " + " | ".join([cell.text for cell in table.rows[0].cells]) + " |")
            markdown_lines.append("|" + "|".join(["---" for _ in table.rows[0].cells]) + "|")

            for row in table.rows[1:]:
                markdown_lines.append("| " + " | ".join([cell.text for cell in row.cells]) + " |")
            markdown_lines.append("\n")

        # Write to output file
        with open(output_path, 'w', encoding='utf-8') as f:
            f.writelines(markdown_lines)

        print(f"✅ Conversion complete! Output: {output_path}")
        return True

    except Exception as e:
        print(f"❌ Error during conversion: {e}")
        import traceback
        traceback.print_exc()
        return False

def extract_images(input_path, output_dir):
    """
    Extract images from Word document

    Args:
        input_path: Path to input Word document
        output_dir: Directory to save extracted images
    """
    try:
        from docx import Document
        import zipfile

        # Create output directory
        os.makedirs(output_dir, exist_ok=True)

        # Word documents are ZIP files containing images
        with zipfile.ZipFile(input_path, 'r') as zip_ref:
            image_files = [f for f in zip_ref.namelist() if f.startswith('media/')]

            if image_files:
                print(f"📁 Extracting {len(image_files)} images to {output_dir}")
                for img_file in image_files:
                    # Extract image
                    img_data = zip_ref.read(img_file)
                    img_name = os.path.basename(img_file)
                    img_path = os.path.join(output_dir, img_name)

                    with open(img_path, 'wb') as f:
                        f.write(img_data)

                print(f"✅ Images extracted successfully!")
            else:
                print("ℹ️  No images found in document")

        return True

    except Exception as e:
        print(f"⚠️  Warning: Could not extract images: {e}")
        return False

def main():
    """Main conversion function"""
    print("=== Word to Markdown Converter (Python Version) ===")
    print("")

    # Check dependencies
    if not check_dependencies():
        sys.exit(1)

    # Define file paths
    input_file = "docs/GLINT tutorial.docx"
    output_file = "docs/tutorial.md"
    images_dir = "docs/images"

    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"❌ Error: Input file not found: {input_file}")
        sys.exit(1)

    # Extract images
    extract_images(input_file, images_dir)

    # Convert to Markdown
    if convert_word_to_markdown(input_file, output_file):
        print("")
        print("=== Conversion Successful! ===")
        print(f"📄 Markdown file: {output_file}")
        print(f"📁 Images directory: {images_dir}")
        print("")
        print("Next steps:")
        print("1. Review the converted Markdown file")
        print("2. Manually fix any formatting issues")
        print("3. Update image references in the Markdown")
        print("4. Test the documentation site: cd docs && jekyll serve")
    else:
        print("")
        print("❌ Conversion failed. Please check the error messages above.")
        sys.exit(1)

if __name__ == "__main__":
    main()
