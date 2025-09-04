# Python 베이스 이미지
FROM python:3.11-slim

# 작업 디렉토리
WORKDIR /app

# 패키지 설치
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 소스 복사
COPY . .

# FastAPI 실행 명령
CMD ["uvicorn", "chemag:app", "--host", "0.0.0.0", "--port", "8000"]
