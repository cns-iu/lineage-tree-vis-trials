# Template for projects using map4sci

## Introduction

The map4sci project template is a [GitHub Template repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-template-repository) that provides a basic template for running map4sci with your network data. Developers can create a new repository based on this template by clicking [here](https://github.com/cns-iu/map4sci-project-template/generate) and using it as a starting point for their own project. The template includes all the necessary files and configurations to run map4sci and includes a run.sh script that can be used to start the map4sci application.

With the map4sci project template, integrating map4sci into your project is easy and straightforward, allowing you to focus on your application logic rather than the infrastructure and configuration required to run map4sci.

## Prerequisites

1. Docker (and Docker Compose) installed and running
2. If running through CWL, cwl-runner must also be installed. You can install the reference runner in python with `python -m pip install cwl-runner cwltool`

## 1. Creating Datasets

To create a dataset, follow these steps:

1. Create a folder with the name of the dataset inside `datasets`.
2. Inside the folder, create a `config.sh` file and a `network.dot` file.
3. For configuring the dataset refer [configuring the dataset](https://github.com/cns-iu/map4sci/tree/code-cleanup-documentation#configuring-the-dataset).

`Note that you can create multiple datasets to be processed.`

## 2. Running

1. `quick run`: Run the `run.sh` script which will loop through all datasets in the datasets directory and run map4sci over them. It skips datasets which have succesfully finished already.
2. `Running through CWL`: Refer [Running using cwltool/cwl-runner](https://github.com/cns-iu/map4sci/tree/code-cleanup-documentation#running-using-cwltoolcwl-runner).
3. `Running through Docker`: Refer [Running using Docker](https://github.com/cns-iu/map4sci/tree/code-cleanup-documentation#running-using-docker).
4. `Running through Docker Compose`: Refer [Running using Docker Compose](https://github.com/cns-iu/map4sci/tree/code-cleanup-documentation#running-using-docker-compose).

## 3. Customization (Advanced)

The script includes several environment variables that can be used to customize the behavior of map4sci. These variables can be set in a constants.sh file located in the same directory as the run.sh script.

To customize the behavior of map4sci, you can use the [constants.sh](https://github.com/cns-iu/map4sci-project-template/blob/main/constants.sh) file and define the variables with the desired values. For example, to specify a different version of map4sci, you can add the following line to the constants.sh file:
`export CURRENT_VERSION=v2`, this changes the version from v1 to v2.

Similarly, you can set the `RUNNER` variable to specify whether to use Docker or CWL as the runner, and set the `CWL_OPTS` variable to specify any additional options to pass to the CWL runner.

If you choose to use CWL, you can also customize the behavior of map4sci by modifying the [map4sci.cwl](https://github.com/cns-iu/map4sci/blob/main/map4sci.cwl) file, which describes the steps and inputs required to generate the map4sci site. You can also use Singularity as a customization option by adding the --singularity flag to the CWL_OPTS variable.

### 3.1 Adding more scripts

In addition to customizing map4sci using environment variables and the constants.sh file, you can also add custom scripts to modify the behavior of map4sci. These scripts should be placed in the scripts directory and should be named with a number prefix to indicate the order in which they should be executed.  For example, a script named 01-my-script.sh will be executed before a script named 02-another-script.sh

To add a custom script, simply create a new file in the scripts directory with the desired name and add the necessary commands to modify the behavior of map4sci. For example, you could add a script to download additional datasets or modify the configuration of the site generator.

Note that the run.sh script will automatically execute any scripts in the scripts directory, in numerical order based on the prefix. This provides a convenient way to add customizations to map4sci without having to modify the core code.

`Note: scripts in the scripts folder should follow a numerical order from 00-*.sh to 99-*.sh (scripts will be skipped if its not of the form NN-*.sh)`

## 4. Results

Results will be in the `raw-data` directory and a fully built static website in the `site` directory.
