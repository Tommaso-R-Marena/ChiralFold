FROM python:3.11-slim-bookworm

RUN apt-get update && apt-get install -y --no-install-recommends \
    libxrender1 libxext6 git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r requirements.txt

COPY app.py .
COPY theme.css .
COPY examples/ examples/

ENV CHIRALFOLD_WEB_HOST=0.0.0.0
ENV CHIRALFOLD_WEB_PORT=7860
EXPOSE 7860

CMD ["python", "app.py"]
