# Template_project



## Getting started

This template is meant as a starting point for creating new repositories for projects in the Industrial Systems Biotechnology group at the institute of Applied Microbiology. This README file will highlight some of the core principles for sharing and creating software/data according to the [FAIR software principles](https://www.nature.com/articles/s41597-022-01710-x).

## A good repository structure
This repository is structured to make it easier to track your work. You can create subdirectories as you please. Make sure you write your code such that it is runnable from the command line from the main directory. This will ensure you do not have any import conflicts.

Please check out what goes where:

**Directories**
- `Data`: Only data as you obtained it. Should never contain anything which is the result of your computations (so also not 'cleaned' data)
- `Figures`: This is where you store the figures which are included in your manuscript AND the scripts to generate them.
- `Models`: All the models used in your work.
- `Results`: Here you can store anything you generated. You can opt to separate several phases of the research (preprocessing, processing, analysis) in several subdirectories.
- `Scripts`: This is where all your notebooks and scripts which you use for a SINGLE analysis. If you use it more often, it is better to store it in the `src` directory.
- `src`: In the `src` directory you can store all your software/reusable code.
- `docs`: For the documentation. the API reference files can be stored in a subdirectory `api_reference`
- `tests`: Here you can store your software tests
- `.github`: This is were github CI action scripts are stored. Please delete if not necessary.

**Files**
- `.gitignore`: every file listed in the `.gitignore` file will be ignored by git. This is very confeniant for automatically generated files and confidential information
- `requirements.txt`: A list of all the packages (including version!) which are required to redo the analyses.
- `.readthedocs.yml`: An example for deploying documentation using readthedocs. Please delete if not necessary.


## Step-by-step plan for using version control with Git
The following step-by-step guidelines will help you create and maintain organized projects on GitLab

1. [Fork this project](https://docs.gitlab.com/ee/user/project/repository/forking_workflow.html). Ensure you create a fork from the main branch and make it private
2. [Clone the project](https://docs.gitlab.com/ee/user/project/repository/index.html#clone-a-repository) to your computer
3. [Maintain the project](https://docs.gitlab.com/ee/user/get_started/get_started_managing_code.html) by doing one commit per change, and pushing regularly to Git. [This webpage](https://carpentries-incubator.github.io/fair-research-software/04-version-control.html) from the software carpentries gives a good overview of how to use version control.
4. If you are working of software development, the following branching structure is recommended:
    - `main`: this is the latest *stable* version of the software
    - `develop`: this is were development happens, only merged into main if it is stable
    - `feature`: in a feature branch, you can develop new features for your software. If the feature is ready, it can be merged into main and the branch can be deleted.

More information about collaborating on software projects can be found [here](https://coderefinery.github.io/git-collaborative/)

## Step-by-step plan for creating repositories using the FAIR principles
The following guidelines will help you build your project such that it is easily Findable, Accessible, Interoperable and Reproducible

1. **Before you start**:
    - Create a venv/conda venv for this specific project
    - Determine which version of Python/pandas/numpy/other basics you will be using and make sure you version pin (e.g. `pandas==0.25.3`) them in your requirements file
    - Replace this [README](https://github.com/hackergrrl/art-of-readme#bonus-the-readme-checklist) with a landing page describing your project and how people should collaborate.

2. **During the development**:
    - Regularly update your requirements file. Preferably after each new installation/deinstallation. Keep in mind that package co-dependencies should be deleted manually to keep a clean list or requirements. with conda, you can use `conda list --from-history` to get all the packages you have installed manually (thus not the co-dependencies). More information and an example can be found in the [code refinery](https://coderefinery.github.io/reproducible-research/dependencies/).
    - See the [Git usage guide](#step-by-step-plan-for-using-version-control-with-git) for how to maintain your repository.
    - Create a test suite in the `tests` directory. Some guidelines can be found from the [carpentries-incubator](https://carpentries-incubator.github.io/fair-research-software/08-code-correctness.html) and [SureSoft](https://suresoft.dev/knowledge-hub/software-testing/). A minimal requirement is writing unit tests using the [Arrange, Act, Assert framework](https://automationpanda.com/2020/07/07/arrange-act-assert-a-pattern-for-writing-good-tests/)
    - Best practice is to introduce a continuous integration (CI) framework. A template can be found in the `.github` folder for [GitHub](https://carpentries-incubator.github.io/fair-research-software/ci-for-testing.html). GitLab works with a `.gitlab-ci.yml` file in the root of the repository. In the [GitLab documentation](https://docs.gitlab.com/ee/tutorials/setup_steps/) and the [SureSoft](https://suresoft.dev/knowledge-hub/continuous-integration/gitlab-ci/) you can find more information about how to set this up. If you want to do something elaborate, you can check the [yaml 'dictionary' from GitLab](https://docs.gitlab.com/ee/ci/yaml/).
    - [Keep writing documentation](https://carpentries-incubator.github.io/better-research-software/07-software-documentation.html)! Even though it is tempting to only do this at the end, it is better to start writing your general documentation describing the principles underlying the software early on. Some good resources can be found at the [GoodDocsProject](https://www.thegooddocsproject.dev/) and the [Diataxis framework](https://diataxis.fr/)
    - Add all confidential files or files/directories such as __pycache__ to the `.gitignore` file. These files will then be ignored by git, and not uploaded to the git server.

3. **Before publication**:
    - **IMPORTANT**: determine the software licence. This tells the user what it is allowed to do with it. Please note that the licence can be affected by potential future patents, or by the packages you use in your software/repository. The EU provides a great tool for finding out which licence is good for you: a [licence compatibility checker](https://interoperable-europe.ec.europa.eu/collection/eupl/solution/licensing-assistant/compatibility-checker) and an [overview of the commonly used licences](https://interoperable-europe.ec.europa.eu/collection/eupl/solution/licensing-assistant/find-and-compare-software-licenses). More information can be found in the [code refinery](https://coderefinery.github.io/social-coding/software-licensing/)
    - Finish your documentation. Make sure you update your [software level documentation](hhttps://carpentries-incubator.github.io/better-research-software/07-software-documentation.html#software-level-documentation) and your [project level documentation](https://carpentries-incubator.github.io/better-research-software/07-software-documentation.html#project-level-documentation)
    - Deploy your documentation using [MKDocs](https://carpentries-incubator.github.io/fair-research-software/09-code-documentation.html#documentation-tools), [ReadTheDocs](https://about.readthedocs.com/?ref=readthedocs.org) (example is in the docs directory), or whatever other tool you'd like to use.
    - Make your repository citable by adding a CITATION file. [This tool](https://citation-file-format.github.io/cff-initializer-javascript/#/) will generate the file for you. More information can be found in the [code refinery](https://coderefinery.github.io/social-coding/software-citation/)
    - You can opt to [publish (versions of) your repository on Zenodo](https://coderefinery.github.io/github-without-command-line/doi/) to get a doi. This would be the FAIRest way to make your software citable.
    - Do you expect people to collaborate or contribute to the repository? Then it is recommended to create a `CONTRIBUTION` file.

4. **Publication**:
    - Create a clean copy of the project on the iAMB GitHub page. 
    - If you do not need the history of your project anymore, [archive it on GitLab](https://docs.gitlab.com/user/project/working_with_projects/#archive-a-project)
    - If you are not working on a confidential project and you do not mind people seeing your git history, you can also [migrate the GitLab repo to GitHub](https://gist.github.com/sxflynn/3ed8f78fe9c4a115ab14857854ab7f6d) without loosing anything.

5. **After publication**: 
    - **Keep on the maintenance!**
    - If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely.
    - Are you often releasing new versions? Then you might want to consider adding a `CHANGELOG` file. The creation of such a file can be automated as described in [this blogpost from Free Code Camp](https://www.freecodecamp.org/news/a-beginners-guide-to-git-what-is-a-changelog-and-how-to-generate-it/).