#!/usr/bin/env python3
"""
Clean up Markdown files converted from Word documents
"""

import re
import sys

def clean_markdown_file(input_file):
    """
    Clean up Markdown file by removing Pandoc-specific formatting

    Args:
        input_file: Path to Markdown file to clean
    """
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            content = f.read()

        # Remove Pandoc image attributes
        content = re.sub(r'\{width="[^"]*"[^}]*\}', '', content)

        # Fix image paths (remove ./docs/ prefix)
        content = re.sub(r'\.\/docs\/images\/', '/GLINT/images/', content)

        # Clean up multiple empty lines
        content = re.sub(r'\n{3,}', '\n\n', content)

        # Fix bold formatting (convert **text** to **text**)
        # Pandoc sometimes creates weird bold patterns
        content = re.sub(r'\*\*\*([^*]+)\*\*\*', r'**\1**', content)

        # Save cleaned content
        with open(input_file, 'w', encoding='utf-8') as f:
            f.write(content)

        print(f"✅ Cleaned: {input_file}")
        return True

    except Exception as e:
        print(f"❌ Error cleaning {input_file}: {e}")
        return False

if __name__ == "__main__":
    if len(sys.argv) > 1:
        file_to_clean = sys.argv[1]
    else:
        file_to_clean = "docs/tutorial.md"

    clean_markdown_file(file_to_clean)
