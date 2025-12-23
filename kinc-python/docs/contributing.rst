Contributing
============

We welcome contributions to KINC-Python! This page explains how to contribute.

Ways to Contribute
------------------

- **Report bugs**: File issues on GitHub
- **Suggest features**: Open feature request issues
- **Improve documentation**: Fix typos, add examples
- **Write code**: Fix bugs, add features
- **Write tests**: Improve test coverage
- **Share use cases**: Tell us how you use KINC-Python

Getting Started
---------------

1. **Fork the repository**:

Visit https://github.com/SystemsGenetics/KINC-Python and click "Fork".

2. **Clone your fork**:

.. code-block:: bash

   git clone https://github.com/YOUR_USERNAME/KINC-Python.git
   cd KINC-Python

3. **Install in development mode**:

.. code-block:: bash

   pip install -e ".[dev]"

This installs KINC-Python plus development tools (pytest, black, flake8).

4. **Create a branch**:

.. code-block:: bash

   git checkout -b feature-name

Development Workflow
--------------------

1. **Make your changes**

2. **Format code**:

.. code-block:: bash

   black kinc/
   flake8 kinc/

3. **Run tests**:

.. code-block:: bash

   pytest tests/

4. **Commit changes**:

.. code-block:: bash

   git add .
   git commit -m "Brief description of changes"

5. **Push to your fork**:

.. code-block:: bash

   git push origin feature-name

6. **Create pull request**:

Visit GitHub and click "New Pull Request".

Code Style
----------

We follow PEP 8 with these conventions:

- **Line length**: 100 characters (not 79)
- **Quotes**: Double quotes for strings
- **Imports**: Grouped (standard, third-party, local)
- **Docstrings**: Google style

Example:

.. code-block:: python

   """
   Module for computing correlations.
   """

   import numpy as np
   from typing import Optional


   def pearson_correlation(x: np.ndarray, y: np.ndarray,
                          min_samples: int = 30) -> float:
       """
       Compute Pearson correlation coefficient.

       Args:
           x: Expression values for first gene
           y: Expression values for second gene
           min_samples: Minimum samples required

       Returns:
           Correlation coefficient or NaN if insufficient data
       """
       # Implementation here
       pass

Use black for automatic formatting:

.. code-block:: bash

   black kinc/

Testing
-------

Write tests for new features:

.. code-block:: python

   # tests/test_correlation.py
   import pytest
   import numpy as np
   from kinc.core.correlation import pearson_correlation


   def test_pearson_perfect_correlation():
       """Test Pearson with perfectly correlated data."""
       x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
       y = np.array([2.0, 4.0, 6.0, 8.0, 10.0])

       corr = pearson_correlation(x, y, min_samples=3)

       assert abs(corr - 1.0) < 1e-10


   def test_pearson_insufficient_samples():
       """Test Pearson with too few samples."""
       x = np.array([1.0, 2.0])
       y = np.array([3.0, 4.0])

       corr = pearson_correlation(x, y, min_samples=5)

       assert np.isnan(corr)

Run tests:

.. code-block:: bash

   pytest tests/
   pytest tests/test_correlation.py -v

Documentation
-------------

Add docstrings to all public functions and classes:

.. code-block:: python

   def function_name(param1: str, param2: int = 10) -> bool:
       """
       Brief one-line description.

       Longer description if needed. Can span multiple lines
       and include additional details about the function.

       Args:
           param1: Description of first parameter
           param2: Description of second parameter (default: 10)

       Returns:
           Description of return value

       Raises:
           ValueError: When param1 is invalid

       Example:
           >>> function_name("test", 5)
           True
       """
       pass

Update RST documentation when adding features:

.. code-block:: bash

   # Edit docs/usage.rst or other files
   cd docs
   make html
   # View in browser: _build/html/index.html

Pull Request Guidelines
------------------------

Good pull requests:

1. **Focused**: Address one issue or feature
2. **Tested**: Include tests for new code
3. **Documented**: Update docstrings and RST docs
4. **Formatted**: Use black and pass flake8
5. **Described**: Clear PR description

PR template:

.. code-block:: text

   ## Description
   Brief description of changes

   ## Type of Change
   - [ ] Bug fix
   - [ ] New feature
   - [ ] Documentation update
   - [ ] Performance improvement

   ## Testing
   How was this tested?

   ## Checklist
   - [ ] Code follows style guidelines
   - [ ] Added tests
   - [ ] Updated documentation
   - [ ] All tests pass

Development Environment
-----------------------

Recommended setup:

.. code-block:: bash

   # Create virtual environment
   python3 -m venv kinc-dev
   source kinc-dev/bin/activate

   # Install development dependencies
   cd KINC-Python
   pip install -e ".[dev]"

   # Install pre-commit hooks (optional)
   pip install pre-commit
   pre-commit install

IDE Setup
~~~~~~~~~

**VS Code**:

.. code-block:: json

   {
     "python.linting.enabled": true,
     "python.linting.flake8Enabled": true,
     "python.formatting.provider": "black",
     "python.formatting.blackArgs": ["--line-length", "100"],
     "editor.formatOnSave": true
   }

**PyCharm**:

- Enable black formatter in preferences
- Set line length to 100
- Enable flake8 linter

Areas Needing Contribution
---------------------------

High Priority
~~~~~~~~~~~~~

- **Parallelization**: MPI, Dask, GPU support
- **Additional analytics**: RMT, powerlaw, corrpower, cond-test
- **Performance optimization**: Cython for hot loops
- **Testing**: Increase test coverage to >80%

Medium Priority
~~~~~~~~~~~~~~~

- **Documentation**: More examples, tutorials
- **Visualization**: Network plotting utilities
- **Benchmarking**: Performance comparisons
- **CI/CD**: Automated testing on push

Low Priority
~~~~~~~~~~~~

- **Notebook examples**: Jupyter notebooks
- **Docker images**: Pre-built containers
- **Web interface**: Optional web GUI

Reporting Issues
----------------

Good bug reports include:

1. **KINC-Python version**: ``kinc --version``
2. **Python version**: ``python --version``
3. **Operating system**: Linux/macOS/Windows
4. **Expected behavior**: What should happen
5. **Actual behavior**: What actually happened
6. **Minimal example**: Code to reproduce
7. **Error message**: Full traceback

Example:

.. code-block:: text

   **Version**: KINC-Python 3.4.2, Python 3.10.0
   **OS**: Ubuntu 22.04

   **Expected**: Similarity should complete without error

   **Actual**: Gets MemoryError after 2 hours

   **Reproduce**:
   ```bash
   kinc similarity --input large.emx.h5 --cmx output.cmx.h5 --ccm output.ccm.h5
   ```

   **Error**:
   ```
   MemoryError: Unable to allocate array...
   ```

Feature Requests
----------------

Good feature requests explain:

1. **Use case**: Why is this needed?
2. **Proposed solution**: How might it work?
3. **Alternatives**: Other approaches considered
4. **Examples**: Concrete usage examples

Code of Conduct
---------------

Be respectful and constructive:

- Use welcoming language
- Accept constructive criticism
- Focus on what's best for the community
- Show empathy towards others

See full `Code of Conduct <https://github.com/SystemsGenetics/KINC-Python/blob/master/CODE_OF_CONDUCT.md>`_.

License
-------

By contributing, you agree that your contributions will be licensed under the same license as KINC-Python.

Questions?
----------

- **Email**: stephen.ficklin@wsu.edu
- **GitHub Discussions**: https://github.com/SystemsGenetics/KINC-Python/discussions
- **Issues**: https://github.com/SystemsGenetics/KINC-Python/issues

Thank You!
----------

Contributors to KINC-Python:

- Stephen Ficklin (maintainer)
- [Your name here - contribute to get listed!]

Special thanks to the original C++ KINC developers:
Ben Shealy, Josh Burns, Frank Feltus, and contributors.
