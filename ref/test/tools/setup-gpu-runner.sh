#!/usr/bin/env bash

# Install drivers
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2404/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt-get update
sudo apt-get -y install cuda-toolkit-12-6
sudo apt-get install nvidia-gds
sudo apt-get install -y cuda-drivers

# reboot

# Install LMOD
sudo apt-get update -y
sudo apt-get install -y lmod
echo "source /usr/share/lmod/lmod/init/bash" >> ~/.bash_profile
source /usr/share/lmod/lmod/init/bash
module list

# Install NVIDIA HPC SDK
curl https://developer.download.nvidia.com/hpc-sdk/ubuntu/DEB-GPG-KEY-NVIDIA-HPC-SDK | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg
echo 'deb [signed-by=/usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' | sudo tee /etc/apt/sources.list.d/nvhpc.list
sudo apt-get update -y
sudo apt-get install -y nvhpc-24-7

# Install cmake
sudo apt-get install -y cmake


# Run persistence driver - not needed?
#sudo /usr/bin/nvidia-persistenced --verbose

# Create a folder
mkdir actions-runner && cd actions-runner
# Download the latest runner package
curl -o actions-runner-linux-x64-2.319.0.tar.gz -L https://github.com/actions/runner/releases/download/v2.319.0/actions-runner-linux-x64-2.319.0.tar.gz
# Optional: Validate the hash
echo "52b8f9c5abb1a47cc506185a1a20ecea19daf0d94bbf4ddde7e617e7be109b14  actions-runner-linux-x64-2.319.0.tar.gz" | shasum -a 256 -c
# Extract the installer
tar xzf ./actions-runner-linux-x64-2.319.0.tar.gz

# Create the runner and start the configuration experience
$ ./config.sh --url https://github.com/NOAA-GSL/SENA-gf --token <given by github interface>
# Last step, run it!
$ ./run.sh
