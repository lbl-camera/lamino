# ReadTheDocs Setup for Tomocam

This directory contains the documentation for Tomocam, built with Sphinx and hosted on ReadTheDocs.

## Building Documentation Locally

### Install Dependencies

```bash
pip install -r requirements.txt
```

### Build HTML Documentation

```bash
cd docs
make html
```

The generated documentation will be in `_build/html/`. Open `_build/html/index.html` in your browser.

### Build PDF Documentation

```bash
make pdf
```

The PDF will be generated in `_build/latex/Tomocam.pdf`.

## Structure

- `index.rst` - Main documentation index
- `getting_started.rst` - Getting started guide
- `installation.rst` - Installation instructions
- `usage.rst` - Usage guide
- `configuration.rst` - Configuration reference
- `examples.rst` - Practical examples
- `contributing.rst` - Contribution guidelines
- `license.rst` - License information
- `conf.py` - Sphinx configuration
- `requirements.txt` - Python dependencies

## ReadTheDocs Configuration

The repository includes `.readthedocs.yaml` at the root, which configures ReadTheDocs to:

1. Install Python dependencies from `docs/requirements.txt`
2. Build the Sphinx documentation
3. Generate PDF, EPUB, and HTML formats

## Automatic Builds

ReadTheDocs will automatically build documentation when:
- Commits are pushed to the main branch
- Pull requests are created
- Tags are created

## Local Development Server

To preview documentation with live reload:

```bash
make serve
```

Then open http://localhost:8000 in your browser.

## Updating Documentation

1. Edit `.rst` files in this directory
2. Rebuild with `make html` to preview
3. Commit and push changes to trigger ReadTheDocs rebuild

## Troubleshooting

**Missing dependencies:**
```bash
pip install -r requirements.txt
```
