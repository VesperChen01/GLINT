#!/usr/bin/env python3
"""
Add Jekyll front matter to Markdown files
"""

import sys

def add_frontmatter(input_file, title="Tutorial"):
    """
    Add Jekyll front matter to a Markdown file

    Args:
        input_file: Path to Markdown file
        title: Page title
    """
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            content = f.read()

        # Check if front matter already exists
        if content.startswith('---'):
            print(f"⚠️  Front matter already exists in {input_file}")
            return True

        # Create front matter
        front_matter = f"""---
layout: default
title: {title}
---

"""

        # Add front matter to the beginning
        new_content = front_matter + content

        # Save the file
        with open(input_file, 'w', encoding='utf-8') as f:
            f.write(new_content)

        print(f"✅ Added front matter to {input_file}")
        return True

    except Exception as e:
        print(f"❌ Error adding front matter to {input_file}: {e}")
        return False

if __name__ == "__main__":
    if len(sys.argv) > 1:
        file_to_process = sys.argv[1]
        title = sys.argv[2] if len(sys.argv) > 2 else "Tutorial"
    else:
        file_to_process = "docs/tutorial.md"
        title = "Tutorial"

    add_frontmatter(file_to_process, title)
