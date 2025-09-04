# 베이스 이미지를 Playwright 공식 이미지로 변경합니다.
# (버전은 필요에 따라 https://mcr.microsoft.com/en-us/product/playwright/python/tags 에서 확인)
FROM mcr.microsoft.com/playwright/python:v1.48.0-jammy

# 작업 디렉토리 설정
WORKDIR /app

# requirements.txt 파일 복사
COPY requirements.txt .

# Playwright는 이미 설치되어 있으므로, 나머지 패키지만 설치합니다.
RUN pip install --no-cache-dir -r requirements.txt

# 소스 코드 복사
COPY . .

# FastAPI 실행 명령
CMD ["uvicorn", "chemag:app", "--host", "0.0.0.0", "--port", "8000"]
