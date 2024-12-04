# Python base project

1. To setup the environnement, make sure to install miniconda https://docs.anaconda.com/miniconda/miniconda-install/
2. Make sur make is installed (for windows : https://gnuwin32.sourceforge.net/packages/make.htm , don't forget to add C:\Program Files (x86)\GnuWin32\bin to Windows' Path)
3. Setup everything by executing the makefile. Open a terminal in this folder and execute :
```
make ENV=my_env
```
Or **for WINDOWS**
```
make -f Makefile_Windows ENV=my_env
```
During this step you will have to enter informations in the terminal. Here is what you have to enter :

* _Please enter the Python interpreter to use :_ Answer 0
* _Project name (python_project_template):_ Write your project name **Avoid spaces and -, it can't be empty**
* _Project version (0.1.0):_ The version you want
* _If yes, it will be installed by default when running `pdm install`. [y/n] (n):_ **Answer y** 
* _Project description (): _ The description you want (can be empty)
* _Which build backend to use?_ **Answer 1 (setuptools)**
* _License(SPDX name) (MIT):_ Leave Empty
* _Author name (Arnaud SOMVILLE):_ Your Name
* _Author email (arnaud.somville@estaca.eu):_ Your email
* _Python requires('*' to allow any) (>=3.11):_ **Leave empty**

## Activating the Environment

To use the project, you need to activate the `my_env` Conda environment created during setup.

### Using VS Code with Conda

1. **Open VS Code:**
   - Navigate to your project folder.

2. **Install the Python Extension:**
   - Open the Extensions view (`Ctrl+Shift+X` on Windows/Linux or `Cmd+Shift+X` on macOS).
   - Search for "Python" and install the official extension by Microsoft.

3. **Select the Environment:**
   - Open the Command Palette (`Ctrl+Shift+P` or `Cmd+Shift+P`).
   - Type `Python: Select Interpreter` and press Enter.
   - Look for the `my_env` environment in the list and select it.

4. **Activate the Environment in the Terminal:**
   - Open the integrated terminal in VS Code (`Ctrl+\` or `Cmd+\`).
   - Run the following command:
     ```bash
     conda activate my_env
     ```

5. **Verify the Environment:**
   - Run the following command in the terminal:
     ```bash
     python --version
     ```
   - Ensure it matches the Python version installed in your `my_env` environment.

## Setting Up Documentation with Sphinx

This project uses Sphinx to generate documentation. Follow the steps below to set up and customize your documentation:

### Steps to Set Up Documentation

1. **Activate the Conda Environment:**
   - Open a terminal in your project directory and activate the `my_env` environment:
     ```bash
     conda activate my_env
     ```

2. **Initialize Sphinx:**
   - Run the Sphinx setup command:
     ```bash
     sphinx-quickstart docs
     ```
   - Follow the prompts to configure Sphinx:
     - **Separate source and build directories?** Yes.
     - **Project name:** Enter your project name.
     - **Author name:** Enter your name.
     - **Project release:** Specify the release version, e.g., `0.1.0`.
     - **Language:** Leave blank for default (English).

3. **Customize `conf.py`:**
   - Open the `docs/conf.py` file to configure Sphinx.
   - Add any required extensions, such as:
     ```python
     extensions = ["sphinx.ext.autodoc", "sphinx.ext.napoleon", "sphinx.ext.viewcode"]
     ```

4. **Generate Documentation:**
   - Build the documentation with:
     ```bash
     make html
     ```
   - The HTML files will be generated in the `docs/_build/html` directory.

5. **Preview the Documentation:**
   - Open the `index.html` file in a web browser:
     ```bash
     open docs/_build/html/index.html
     ```

### Customizing the Documentation

- **Add Content:**
  - Write your documentation in `.rst` files within the `docs/source` directory.
  - Update the `index.rst` file to include references to your new files.

- **Automate API Documentation:**
  - Use `sphinx.ext.autodoc` to automatically generate documentation from your Python code:
    ```bash
    sphinx-apidoc -o docs/source/ src/
    ```

With these steps, you can easily set up and manage your project documentation using Sphinx.

## Deploying Documentation to GitHub Pages

This project is set up to automatically generate and deploy Sphinx documentation to GitHub Pages using a GitHub Actions workflow. Follow these simple steps to enable deployment of your documentation:

### Steps to Enable Documentation Deployment

1. **Ensure GitHub Pages is Enabled:**
   - Go to your repository on GitHub.
   - Navigate to `Settings > Pages`.
   - Under the "Source" section, select the `gh-pages` branch as the source for your GitHub Pages. 
   - Leave the folder as `/ (root)`.
   - Click **Save**.

2. **Trigger the Workflow:**
   - Push changes to the `main` branch or create a pull request. The CI will:
     1. Build your Sphinx documentation.
     2. Deploy the HTML files to the `gh-pages` branch.

3. **Access Your Documentation:**
   - Once deployed, your documentation will be accessible at:
     ```
     https://<your-username>.github.io/<your-repository-name>/
     ```
   - Replace `<your-username>` and `<your-repository-name>` with your GitHub username and repository name.

### Troubleshooting

- If the pages do not appear:
  - Verify that the `gh-pages` branch exists and contains your HTML files.
  - Confirm that GitHub Pages is enabled in your repository settings.
  - Check the workflow logs in the Actions tab for errors.
With these steps, your Sphinx documentation will be live and ready to share!

- Environment Activation Fails:
  - Ensure Miniconda is properly installed.
  - Verify the environment exists using conda env list. If it doesn’t, rerun make ENV=my_env.
