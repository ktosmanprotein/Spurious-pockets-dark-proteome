FROM python:3.11-slim
WORKDIR /app
COPY *.py ./
CMD ["python", "audit_numbers_p2.py"]
