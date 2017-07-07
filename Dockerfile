FROM python:2.7

# Install / update relevant ubuntu packages:
RUN apt-get update \
	&& apt-get install -y --no-install-recommends build-essential unzip wget libgmp3-dev python-pip

# Download and install glpk:
RUN mkdir /usr/local/glpk \
	&& curl http://ftp.gnu.org/gnu/glpk/glpk-4.39.tar.gz \
	| tar xvzC /usr/local/glpk --strip-components=1 \
	&& cd /usr/local/glpk \
	&& ./configure \
	&& make \
	&& make install
	
# Install pyglpk:
RUN pip install --upgrade pip \
	&& pip install setuptools \
	&& pip install glpk
	
# Install libSBML
RUN pip install python-libsbml

# Install Cobra toolbox
RUN pip install cobra

#Flask and Flask-RESTful
RUN pip install flask \
	&& pip install flask_restful

# Update paths:
ENV LD_LIBRARY_PATH /usr/local/lib:${LD_LIBRARY_PATH}

# Install subliminal-py
RUN pip install subliminal-py

# Make current directory visible inside Docker container:
RUN cd home/
#COPY rest.py /home
WORKDIR /home
RUN git clone https://github.com/dbkgroup/reaction-balancer.git
RUN ls

# Run test:
ENTRYPOINT ["python"]
CMD ["reaction-balancer/rest.py"]