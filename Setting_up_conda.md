# Creating a conda environment

First, we will need to install Miniforge.
 - Please follow the instructions below which are obtained from miniforge [GitHub repository](https://github.com/conda-forge/miniforge?tab=readme-ov-file#unix-like-platforms-macos-linux--wsl)

Download the installation bash script.
```
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
```

Run the script to install Miniforge.
```
bash Miniforge3-$(uname)-$(uname -m).sh
```

If you are using MacOS, install Miniforge by first installing Homebrew:
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Once already installed, install Miniforge as follows:
```
brew install miniforge
```

Finally run this command:
```
conda init zsh
```
Restart your terminal.

Following successful Miniforge installation, we can use `conda` to install the tools we need to an environment in your local computer.
Install the tools as follows:

```
conda env create -f environment.yml
conda env create -f snippy-environment.yml
conda env create -f artic-mpox-environment.yml
```

### Creating a Python environment
We also need to create a python environment to install Python packages that are needed by some python scripts we will run in some analysis steps.
```
python3 -m venv ./py3env
source ./py3env/bin/activate
python3 -m pip install -r ./artic-requirements.txt
```
