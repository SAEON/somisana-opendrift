FROM opendrift/opendrift:latest

# we just need to add our configuration code to the standard opendrift image 
WORKDIR /somisana
ADD . /somisana

# Install somisana-opendrift
RUN pip install -e .

ENTRYPOINT ["python", "cli.py"]
