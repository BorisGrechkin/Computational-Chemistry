FROM python:3.10.5-slim-buster

RUN mkdir -p /Results/Coordinates
RUN mkdir -p /Results/GIF_files

COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt

COPY Scripts/tools.py Scripts/tools.py
COPY convert.py convert.py

ENTRYPOINT ["python", "convert.py"]
