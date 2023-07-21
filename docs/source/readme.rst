.. image:: https://user-images.githubusercontent.com/25851824/200928507-a65327f9-bc70-4c12-beaa-bc6ba74d968e.svg
   :alt: logo
   :width: 60%
   :align: center


`Docs <https://bessagroup.github.io/bessa-pypi-template/>`__ | `GitHub <https://github.com/bessagroup/bessa-pypi-template>`__

.. <!-- | `Installation <link_to_installation_instructions>` -->
.. <!-- | `PyPI <link_to_pypi_package_website>` -->

**First publication:** June 14, 2023

Summary
-------

This repository serves as a template for Python code. The template is compliant to the `Bessa Research Group Python Development Code of Conduct <https://github.com/bessagroup/python_code_of_conduct>`__.

The repository is suitable for any Python code that works in version 3.7+.

Statement of need
-----------------

Members of the Bessa Research Group can use this template to create new Python repositories. The template is compliant to the `Bessa Research Group Python Development Code of Conduct <https://github.com/bessagroup/python_code_of_conduct>`__.

Authorship
----------

**Authors**:
    - Martin van der Schelling (M.P.vanderSchelling@tudelft.nl)

**Authors afilliation:**
    - Bessa Research Group @ Delft University of Technology

**Maintainer:**
    - Martin van der Schelling (M.P.vanderSchelling@tudelft.nl)

**Maintainer afilliation:**
    - Bessa Research Group @ Delft University of Technology

Getting started
---------------

Using this as a template for a new repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When you create a new repository through the GitHub website, you can select this template as a starting point. This will create a new repository with the same structure as this one. You can then clone the repository to your local machine and start working on your code and replacing the template code with your own code.

Using this as a template for an existing repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have an existing repository, you can clone the repository to your local machine and copy the necessary files over manually.

GitHub pages functionality
~~~~~~~~~~~~~~~~~~~~~~~~~~

If you want to use the GitHub Pages functionalities of this template repository, please check the box to clone all branches. The sphinx documentation is build automatically in the ``gh-pages`` branch with GitHub Workflows whenever you push to the ``main`` branch.

In order to commit and push to the ``gh-pages`` branch, you need to authenticate the workflow with your GitHub credentials. This is done with GitHub secrets. This is not done automatically when you clone the repository, so you need to do this manually:

1. Go to your GitHub repository page and click on the "Settings" tab.
2. Click on "Secrets" in the left sidebar menu.
3. Click on the "New secret" button.
4. Enter ``GHPAGES_TOKEN`` as the name of the secret.
5. Generate a new token by clicking on the "Generate a new token" link.
6. Give the token a name and select the appropriate scopes.
7. Click on the "Generate token" button.
8. Copy the generated token and paste it into the "Value" field.
9. Click on the "Add secret" button to save the token.

Once you've created the secret ``GHPAGES_TOKEN``, it can be used in your GitHub workflow scripts by referencing it using the ``${{ secrets.GHPAGES_TOKEN }}`` syntax. This provides a secure way to authenticate with the GitHub API and perform actions such as pushing to a repository, creating issues, and deploying to GitHub Pages.

More about GitHub secrets can be found `here <https://docs.github.com/en/actions/security-guides/encrypted-secrets>`__.

Community Support
-----------------

If you find any issues, bugs or problems with this template, please use the `GitHub issue tracker <https://github.com/bessagroup/bessa-pypi-template/issues>`__ to report them.

License
-------

Copyright 2023, Martin van der Schelling

All rights reserved.

This project is licensed under the BSD 3-Clause License. See `LICENSE <https://github.com/bessagroup/bessa-pypi-template/blob/main/LICENSE>`__ for the full license text.
