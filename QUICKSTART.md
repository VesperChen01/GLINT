# 🚀 GLINT Documentation Quick Start Guide

This guide helps you quickly set up professional GitHub Pages documentation for GLINT using your Word tutorial.

## 📋 What You'll Get

✅ **Professional documentation website** with GitHub Pages
✅ **Converted tutorial** from Word to Markdown
✅ **Modern, responsive design** with custom styling
✅ **Auto-deployment** via GitHub Actions
✅ **Complete documentation structure** (Home, Installation, Tutorial, API)

## 🎯 Two Simple Steps

### Step 1: Convert Your Word Document

Choose **ONE** of the following methods:

#### Method A: Using Pandoc (Recommended - Best Quality)

```bash
# Install pandoc (one-time setup)
brew install pandoc  # macOS
# or
sudo apt-get install pandoc  # Linux
# or download from https://pandoc.org/installing.html (Windows)

# Run conversion
bash convert_docx.sh
```

#### Method B: Using Python (Alternative)

```bash
# Install python-docx (one-time setup)
pip install python-docx

# Run conversion
python3 convert_docx_python.py
```

### Step 2: Set Up Documentation Website

```bash
# Run the setup script
bash setup_github_pages.sh

# Test locally (optional)
cd docs
jekyll serve
```

## 📁 What Gets Created

```
GLINT/
├── docs/
│   ├── _config.yml          # Jekyll configuration
│   ├── _layouts/
│   │   └── default.html     # Page layout
│   ├── assets/
│   │   └── css/
│   │       └── custom.css   # Custom styling
│   ├── images/              # Extracted images
│   ├── index.md             # Home page
│   ├── installation.md      # Installation guide
│   ├── tutorial.md          # Your converted tutorial
│   ├── api.md               # API reference
│   └── .nojekyll            # GitHub Pages file
├── .github/
│   └── workflows/
│       └── deploy-docs.yml  # Auto-deployment
├── convert_docx.sh          # Pandoc converter
├── convert_docx_python.py   # Python converter
└── setup_github_pages.sh    # Setup script
```

## 🌐 Deploy to GitHub Pages

### 1. Commit Your Changes

```bash
git add docs/ .github/workflows/ *.sh *.py
git commit -m "Add GitHub Pages documentation"
git push
```

### 2. Enable GitHub Pages

1. Go to your repository on GitHub
2. Navigate to **Settings** > **Pages**
3. Under **Source**, select **GitHub Actions**
4. Save changes

### 3. Watch Deployment

1. Go to the **Actions** tab
2. Click on "Deploy Documentation" workflow
3. Wait for the green checkmark ✅

### 4. Visit Your Site

Your documentation will be live at:
```
https://<your-username>.github.io/GLINT/
```

## 🎨 Customize Your Documentation

### Edit Content

Open and edit these files:
- `docs/index.md` - Home page content
- `docs/tutorial.md` - Your tutorial (after conversion)
- `docs/installation.md` - Installation instructions
- `docs/api.md` - API documentation

### Change Styling

Edit `docs/assets/css/custom.css`:
```css
:root {
  --primary-color: #3b82f6;  /* Change main color */
  --secondary-color: #10b981; /* Change accent color */
}
```

### Update Navigation

Edit `docs/_config.yml`:
```yaml
nav:
  - title: Home
    url: /
  - title: Installation
    url: /installation
  - title: Tutorial
    url: /tutorial
  - title: API Reference
    url: /api
```

## 🔧 Troubleshooting

### Conversion Issues

**Problem:** Word document doesn't convert well
**Solution:**
1. Try both conversion methods (pandoc vs python)
2. Manually edit the resulting Markdown file
3. Copy-paste content from Word to Markdown directly

**Problem:** Images not extracted
**Solution:**
1. Manually save images from Word to `docs/images/`
2. Update image references in Markdown: `![alt](/GLINT/images/filename.png)`

### Deployment Issues

**Problem:** GitHub Actions fails
**Solution:**
1. Check repository Settings > Pages > Source is "GitHub Actions"
2. Verify workflow file syntax is correct
3. Check Actions tab for specific error messages

**Problem:** Site shows 404 error
**Solution:**
1. Wait a few minutes for deployment to complete
2. Check the URL is correct: `https://<username>.github.io/GLINT/`
3. Verify `.nojekyll` file exists in `docs/`

## 📊 Comparison of Conversion Methods

| Feature | Pandoc | Python (docx) |
|---------|--------|---------------|
| **Formatting Quality** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ |
| **Image Extraction** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ |
| **Table Support** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ |
| **Ease of Use** | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| **Installation** | Moderate | Easy |

**Recommendation:** Try Pandoc first for best quality, fall back to Python if needed.

## 🎓 Learning Resources

- [Markdown Guide](https://www.markdownguide.org/) - Learn Markdown syntax
- [Jekyll Documentation](https://jekyllrb.com/docs/) - Understand Jekyll
- [GitHub Pages Guide](https://docs.github.com/en/pages) - GitHub Pages documentation
- [Bootstrap 5](https://getbootstrap.com/docs/5.1/) - Styling framework used

## 💡 Pro Tips

1. **Test Locally First:** Always run `jekyll serve` to test changes before pushing
2. **Keep Backups:** Commit frequently to avoid losing work
3. **Use Relative Paths:** For images and links, use `/GLINT/path/to/file`
4. **Mobile Friendly:** Test your documentation on mobile devices
5. **SEO Friendly:** The setup includes SEO tags for better search visibility

## 🆘 Need Help?

- **Documentation Issues:** Check [DOCUMENTATION_SETUP.md](DOCUMENTATION_SETUP.md)
- **GLINT Issues:** [GitHub Issues](https://github.com/VesperChen01/GLINT/issues)
- **General Questions:** [GitHub Discussions](https://github.com/VesperChen01/GLINT/discussions)

## ✅ Checklist

Before deploying, make sure:

- [ ] Word document converted to Markdown
- [ ] Images extracted and referenced correctly
- [ ] All links work in the Markdown
- [ ] Tested locally with `jekyll serve`
- [ ] Committed all changes to git
- [ ] GitHub Pages enabled in repository settings
- [ ] GitHub Actions workflow is active

## 🎉 You're Ready!

Once deployed, your professional documentation will help users:
- Understand GLINT features quickly
- Follow step-by-step tutorials
- Reference API documentation
- Get support easily

Your documentation website will automatically update when you push changes to the `docs/` directory!

---

**Happy documenting! 📚✨**
