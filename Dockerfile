FROM ghcr.io/saeon/opendrift_v1.11.0

# we just need to add our configuration code to the standard opendrift image 
WORKDIR /somisana
ADD . /somisana

# Install somisana-opendrift
RUN pip install -e .

# run cli.py as the default behaviour
ENTRYPOINT ["python", "cli.py"]
