FROM continuumio/miniconda3
SHELL ["/bin/bash", "-c"]

RUN apt-get update
RUN apt-get install -y libfontconfig1 libxrender1

# get and install vina
RUN mkdir vina
WORKDIR /vina
RUN wget http://vina.scripps.edu/download/autodock_vina_1_1_2_linux_x86.tgz
RUN tar xzvf autodock_vina_1_1_2_linux_x86.tgz

# Python packages from conda environment.yml file
ADD environment.yml /tmp/environment.yml
RUN chmod -R 777 /tmp

# add conda to bashrc
RUN echo 'export PATH="/opt/conda/bin:$PATH"' >> ~/.bashrc
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc

# source bashrc
RUN source ~/.bashrc

# create conda env from environment.yml
ADD environment.yml /tmp/environment.yml
WORKDIR /tmp
RUN conda env create -f environment.yml
WORKDIR /
RUN rm -rf /tmp