#!/bin/bash

sudo docker build -t ghcr.io/ahmedmoustafa/bioinformatics-tools .
# sudo docker run -v /data:/data -it ghcr.io/ahmedmoustafa/bioinformatics-tools /bin/bash
# ~/token.txt is a Github personal access tokens
cat ~/token.txt | sudo docker login ghcr.io -u ahmedmoustafa --password-stdin
sudo docker push ghcr.io/ahmedmoustafa/bioinformatics-tools
