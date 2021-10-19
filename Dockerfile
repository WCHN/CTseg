FROM ubuntu:20.04

MAINTAINER Guillaume Flandin <g.flandin@ucl.ac.uk>

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get -y install \
     unzip xorg wget \
 && apt-get clean \
 && rm -rf \
     /tmp/hsperfdata* \
     /var/*/apt/*/partial \
     /var/lib/apt/lists/* \
     /var/log/apt/term*

# Install MATLAB MCR in /opt/mcr/
RUN mkdir /opt/mcr_install
RUN mkdir /opt/mcr
RUN wget --progress=bar:force -P /opt/mcr_install https://ssd.mathworks.com/supportfiles/downloads/R2021b/Release/0/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_R2021b_glnxa64.zip
RUN unzip -q /opt/mcr_install/MATLAB_Runtime_R2021b_glnxa64.zip -d /opt/mcr_install
RUN /opt/mcr_install/install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent
RUN rm -rf /opt/mcr_install /tmp/*

# Install SPM Standalone in /opt/spm12/
ENV LD_LIBRARY_PATH /opt/mcr/v911/runtime/glnxa64:/opt/mcr/v911/bin/glnxa64:/opt/mcr/v911/sys/os/glnxa64:/opt/mcr/v911/sys/opengl/lib/glnxa64:/opt/mcr/v911/extern/bin/glnxa64
ENV MCR_INHIBIT_CTF_LOCK 1
ENV SPM_HTML_BROWSER 0
# Running SPM once with "function exit" tests the succesfull installation *and*
# extracts the ctf archive which is necessary if singularity is going to be
# used later on, because singularity containers are read-only.
# Also, set +x on the entrypoint for non-root container invocations
RUN wget --no-check-certificate --progress=bar:force -P /opt https://www.fil.ion.ucl.ac.uk/spm/download/restricted/utopia/dev/soon_gone/spm12_r8168_BI_Linux_R2021b.zip
RUN unzip -q /opt/spm12_r8168_BI_Linux_R2021b.zip -d /opt
RUN rm -f /opt/spm12_r8168_BI_Linux_R2021b.zip
RUN /opt/spm12/spm12 function exit
RUN chmod +x /opt/spm12/spm12

# Hack to ensure the CTseg model files are only downloaded once
RUN /opt/spm12/spm12 eval "spm_CTseg(1)"; exit 0

# Configure entry point
ENTRYPOINT ["/opt/spm12/spm12"]
CMD ["--help"]