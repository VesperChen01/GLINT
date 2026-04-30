#!/bin/bash
# Word to Markdown Converter Script
# This script helps convert Word documents to Markdown format

echo "=== Word to Markdown Converter ==="
echo ""

# Check if pandoc is installed
if command -v pandoc &> /dev/null; then
    echo "✅ Pandoc found! Using pandoc for conversion..."
    echo ""

    # Convert the tutorial document
    INPUT_FILE="docs/GLINT tutorial.docx"
    OUTPUT_FILE="docs/tutorial.md"

    if [ -f "$INPUT_FILE" ]; then
        echo "Converting: $INPUT_FILE"
        pandoc "$INPUT_FILE" -o "$OUTPUT_FILE" --extract-media=./docs/images/
        echo "✅ Conversion complete! Output: $OUTPUT_FILE"
        echo "📁 Images extracted to: docs/images/"
    else
        echo "❌ Error: Input file not found: $INPUT_FILE"
        exit 1
    fi

else
    echo "❌ Pandoc not found. Please install it first:"
    echo ""
    echo "macOS: brew install pandoc"
    echo "Ubuntu/Debian: sudo apt-get install pandoc"
    echo "Windows: Download from https://pandoc.org/installing.html"
    echo ""
    echo "Or install python-docx:"
    echo "pip install python-docx"
    echo ""
    exit 1
fi

echo ""
echo "=== Next Steps ==="
echo "1. Review the converted Markdown file: docs/tutorial.md"
echo "2. Clean up the formatting if needed"
echo "3. Run setup_github_pages.sh to create the documentation website"
