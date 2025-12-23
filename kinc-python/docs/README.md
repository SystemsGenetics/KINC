# KINC-Python Documentation

This directory contains Sphinx documentation for KINC-Python.

## Building Documentation

### Install Dependencies

```bash
pip install -r requirements.txt
```

### Build HTML Documentation

```bash
make html
```

The built documentation will be in `_build/html/index.html`.

### Build PDF Documentation

```bash
make latexpdf
```

### Clean Build

```bash
make clean
```

## Documentation Structure

- `index.rst` - Main documentation index
- `about.rst` - Project overview and background
- `installation.rst` - Installation instructions
- `quickstart.rst` - 5-minute quick start guide
- `usage.rst` - Detailed command-line usage
- `api.rst` - Python API reference
- `examples.rst` - Real-world examples and workflows
- `data_formats.rst` - HDF5 file format documentation
- `performance.rst` - Performance tuning guide
- `troubleshooting.rst` - Common issues and solutions
- `cite.rst` - Citation information
- `contributing.rst` - Contribution guidelines

## Publishing to ReadTheDocs

1. Create account at https://readthedocs.org
2. Import KINC-Python repository
3. Build will automatically trigger on commits
4. Documentation will be available at https://kinc-python.readthedocs.io

## Local Preview

```bash
# Build and serve locally
make html
python3 -m http.server 8000 --directory _build/html

# Visit http://localhost:8000 in browser
```

## Style Guide

- Use reStructuredText (RST) format
- Follow Sphinx conventions
- Include code examples
- Cross-reference other sections with `:doc:`
- Use proper heading hierarchy (=, -, ~, ^)
