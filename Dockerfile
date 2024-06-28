FROM opendrift/opendrift:latest

# we just need to add our configuration code to the standard opendrift image 
WORKDIR /somisana
ADD . /somisana

# Install somisana-opendrift
RUN pip install -e .

# run cli.py as the default behaviour
# (I had to use the full path to the python exe otherwise it wouldn't find it, which I don't understand as /opt/conda/bin is in $PATH?)
ENTRYPOINT ["/opt/conda/bin/python", "cli.py"]
