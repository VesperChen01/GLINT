# GLINT Documentation Setup Guide

This guide will help you set up a professional GitHub Pages documentation website for GLINT using your Word tutorial document.

## Prerequisites

### Option 1: Using Pandoc (Recommended)

Pandoc is the most reliable tool for converting Word documents to Markdown.

**macOS:**
```bash
brew install pandoc
```

**Ubuntu/Debian:**
```bash
sudo apt-get install pandoc
```

**Windows:**
Download from [https://pandoc.org/installing.html](https://pandoc.org/installing.html)

### Option 2: Using Python

If you prefer Python tools:

```bash
pip install python-docx
```

## Step 1: Convert Word Document to Markdown

Run the conversion script:

```bash
bash convert_docx.sh
```

This will:
1. Convert `docs/GLINT tutorial.docx` to `docs/tutorial.md`
2. Extract images to `docs/images/`
3. Preserve formatting and structure

### Manual Conversion (if script fails)

If the automated conversion doesn't work well, you can:

1. **Open your Word document**
2. **Copy the content** section by section
3. **Paste into** `docs/tutorial.md`
4. **Convert formatting:**
   - Headings: Use `#`, `##`, `###` for H1, H2, H3
   - Bold: `**text**`
   - Italic: `*text*`
   - Code: \`code\` or ```code```
   - Links: `[text](url)`
   - Images: `![alt](path/to/image)`

## Step 2: Set Up GitHub Pages Infrastructure

Run the setup script:

```bash
bash setup_github_pages.sh
```

This creates:
- ✅ Jekyll configuration (`docs/_config.yml`)
- ✅ Custom styling (`docs/assets/css/custom.css`)
- ✅ Page layouts (`docs/_layouts/default.html`)
- ✅ Documentation pages (Home, Installation, Tutorial, API)
- ✅ GitHub Actions workflow for auto-deployment
- ✅ `.nojekyll` file for proper GitHub Pages handling

## Step 3: Test Locally

### Install Jekyll and dependencies

**macOS/Linux:**
```bash
# Install Ruby (if not already installed)
# macOS: brew install ruby
# Ubuntu: sudo apt-get install ruby-full build-essential zlib1g-dev

# Install Jekyll and gems
gem install jekyll jekyll-seo-tag jekyll-sitemap
```

**Windows:**
Follow the [Jekyll on Windows](https://jekyllrb.com/docs/installation/windows/) guide.

### Preview your site

```bash
cd docs
jekyll serve
```

Your site will be available at: `http://localhost:4000/GLINT/`

## Step 4: Customize Your Content

### Edit the tutorial

Open `docs/tutorial.md` and:
1. Review the converted content
2. Fix any formatting issues
3. Add screenshots and diagrams
4. Ensure all links work correctly

### Customize styling

Edit `docs/assets/css/custom.css` to:
- Change color scheme
- Adjust fonts and spacing
- Add custom components
- Modify responsive behavior

### Update configuration

Edit `docs/_config.yml` to:
- Change site title and description
- Update navigation menu
- Add plugins
- Configure theme settings

## Step 5: Deploy to GitHub Pages

### Enable GitHub Pages

1. Go to your repository on GitHub
2. Navigate to **Settings** > **Pages**
3. Under **Source**, select **GitHub Actions**
4. The workflow we created will handle the rest

### Push your changes

```bash
git add docs/ .github/workflows/
git commit -m "Add GitHub Pages documentation"
git push
```

### Monitor deployment

1. Go to **Actions** tab in your repository
2. Watch the "Deploy Documentation" workflow
3. Once complete, your site will be live at:
   `https://<username>.github.io/GLINT/`

## Step 6: Custom Domain (Optional)

If you want to use a custom domain:

1. Buy a domain (e.g., `glint-project.org`)
2. Add a `CNAME` file in `docs/`:
   ```
   glint-project.org
   ```
3. Configure DNS settings with your domain provider
4. Update GitHub Pages settings with your custom domain

## Maintenance

### Updating documentation

```bash
# 1. Make changes to Markdown files
# 2. Test locally
cd docs && jekyll serve

# 3. Commit and push
git add docs/
git commit -m "Update documentation"
git push
```

### Adding new pages

1. Create a new Markdown file in `docs/`
2. Add front matter:
   ```yaml
   ---
   layout: default
   title: Your Page Title
   ---
   ```
3. Add to navigation in `docs/_config.yml`
4. Test and deploy

## Troubleshooting

### Pandoc not found

**Error:** `pandoc: command not found`

**Solution:** Install pandoc using the instructions in Prerequisites section.

### Jekyll build errors

**Error:** Various Ruby/gem errors

**Solution:**
```bash
# Update Ruby gems
gem update --system

# Reinstall Jekyll
gem uninstall jekyll
gem install jekyll
```

### Images not displaying

**Error:** Images broken on GitHub Pages

**Solution:**
1. Ensure images are in `docs/images/`
2. Use correct relative paths: `![alt](/GLINT/images/filename.png)`
3. Check image filenames (case-sensitive)

### GitHub Actions deployment fails

**Error:** Workflow fails in Actions tab

**Solution:**
1. Check the workflow logs for specific errors
2. Ensure repository has GitHub Pages enabled
3. Verify workflow file syntax is correct
4. Check that all required permissions are granted

## Advanced Customization

### Adding a search function

Install `jekyll-search` plugin:

```bash
gem install jekyll-search
```

Add to `docs/_config.yml`:
```yaml
plugins:
  - jekyll-search
```

### Adding analytics

Add Google Analytics or other tracking to `docs/_layouts/default.html`:

```html
<!-- Google Analytics -->
<script async src="https://www.googletagmanager.com/gtag/js?id=GA_MEASUREMENT_ID"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());
  gtag('config', 'GA_MEASUREMENT_ID');
</script>
```

### Custom themes

Replace the default theme with a custom Jekyll theme:

```yaml
# In docs/_config.yml
theme: minima  # or any other Jekyll theme
```

## Resources

- [Jekyll Documentation](https://jekyllrb.com/docs/)
- [GitHub Pages Guide](https://docs.github.com/en/pages)
- [Markdown Guide](https://www.markdownguide.org/)
- [Pandoc Manual](https://pandoc.org/MANUAL.html)

## Support

If you encounter issues:

1. Check the [GitHub Issues](https://github.com/VesperChen01/GLINT/issues)
2. Review Jekyll and GitHub Pages documentation
3. Search for similar problems online
4. Ask for help in GitHub Discussions

## Next Steps

1. ✅ Convert your Word document
2. ✅ Set up the documentation infrastructure
3. ✅ Test locally
4. ✅ Deploy to GitHub Pages
5. ✅ Share your documentation with users!

Your professional documentation website will help users understand and use GLINT effectively!
