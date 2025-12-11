## Development Setup

To set up a development environment for this project, follow these steps:

1. **Create a new conda environment** (replace `cnit` with your preferred name if desired):

   ```bash
   conda create --name cnit python=3.9
   ```

2. **Activate your environment**:

   ```bash
   conda activate cnit
   ```

3. **Install all dependencies and set up the package**:

   ```bash
   make env
   ```

This will install all required dependencies (including notebook, docs, and testing tools) and install the package in editable mode.

> **Note:**  
> The Makefile will check that you are not in the `base` environment and will refuse to run if you are.

---

**TL;DR:**
```bash
conda create --name cnit python=3.9
conda activate cnit
make env
```