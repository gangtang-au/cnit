Contributing
============

CNit is developed and maintained by Gang Tang as part of the MAGICC project.

How to Contribute
-----------------

We welcome contributions!  If you're interested in contributing:

1. **Report bugs:** Open an issue on `GitHub <https://github.com/gangtang-au/cnit/issues>`_
2. **Suggest features:** Describe your idea in a GitHub issue
3. **Ask questions:** Email gang.tang.au@gmail.com or open a discussion

For substantial contributions or collaborations, please contact Gang Tang directly
to discuss coordination.

Development Setup
-----------------

If you want to work with CNit's source code:

..  code-block::  bash

   # Clone the repository
   git clone https://github.com/gangtang-au/cnit.git
   cd cnit

   # Create conda environment
   conda create --name cnit python=3.9
   conda activate cnit

   # Install in development mode
   make env

   # Run tests
   make test

   # Build documentation
   make docs

License
-------

CNit is released under the BSD 3-Clause License.  See the LICENSE file for details.