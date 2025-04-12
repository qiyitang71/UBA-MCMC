# Use the official Ubuntu base image
FROM ubuntu:24.04

# Set the maintainer label
LABEL maintainer="qiyi.tang71@gmail.com"


#export timezone
ENV TZ=Europe/London

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y curl
RUN curl -I http://archive.ubuntu.com/ubuntu/
# Update and install basic dependencies
RUN apt-get update && \
    apt-get install -y \
    wget \
    gnupg \
    curl \
    lsb-release \
    ca-certificates \
    software-properties-common \
    python3 \
    python3-pip \
    openjdk-11-jdk \
    && apt-get clean

# Install Spot from official repo
RUN curl -fsSL https://www.lrde.epita.fr/repo/debian.gpg -o /etc/apt/trusted.gpg.d/spot.gpg && \
    echo 'deb http://www.lrde.epita.fr/repo/debian/ stable/' >> /etc/apt/sources.list && \
    apt-get update && \
    apt-get install -y spot libspot-dev && \
    apt-get clean

# Set environment variables for Java
ENV JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
ENV PATH="$JAVA_HOME/bin:$PATH"


# Set the working directory
WORKDIR /workspace

# Copy your source code into the container
# Replace 'src' with your actual source code directory or file
COPY UBA-MCMC/ /workspace/

# Default command (optional)
CMD ["bash"]
