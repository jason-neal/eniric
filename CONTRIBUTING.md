## Contributing to ENIRIC

:+1::tada: First off, thanks for taking the time to contribute! :tada::+1:


#### Do you have more general questions about eniric?

* Please check the [documentation](https://eniric.readthedocs.io/en/latest) to make sure your question is not answered there.
* If not, feel free to open a new issue

#### **Did you find a bug?**

* Please open a GitHub [issue](https://github.com/jason-neal/eniric/issues). to tell us about it.
  Would be even more awesome if you could provide a possible solution!  
  &nbsp; maybe check if the bug was not already reported (and solved) by searching the closed issues.

* If possible, use the bug report template to create the issue.


#### Did you write a patch that fixes a bug?
* Did you run the tests locally/do they pass? You can use `make test` from the root directory.
* Did you enable pre-commit hooks?
* Open a new GitHub pull request.
* Ensure the pull request clearly describes the problem and solution. Include the relevant issue number if applicable.


#### Coding conventions
The style guide for eniric is governed by [black](https://github.com/ambv/black) wth default settings.
You can run this yourself but it is also included in the pre-commit hooks.

Before commiting changes you can activate the [pre-commit hooks](https://github.com/pre-commit/pre-commit) using  

    pip install pre-commit
    pre-commit install

This will then preform the pre-commit checks before commiting.

Thank you!
