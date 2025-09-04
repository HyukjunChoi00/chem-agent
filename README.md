# ChemAgent

- Google Gemini
- FastAPI
- RDkit tool

### Usage Example
- Draw CCCCC
- give me the smiles representation of testosterone and draw it please
- what is the smiles representation of hexane?

```python
pip install -r requirements.txt
```

### Run
```
python chemag.py
```

### Finish
```
taskkill /F /IM python.exe
```

### docker image build
```
git clone https://github.com/HyukjunChoi00/chem-agent
```

```
cd chem-agent
```

```
docker build -t chem-agent-app .
```

위 코드까지 입력 시 docker image를 만들 수 있으며, AWS Elastic Container Registry에 image를 push하여 AWS App Runner와 같은 곳에 **배포**가 가능합니다.  
```
docker run -p 8000:8000 -e GOOGLE_API_KEY="YOUR API KEY" chem-agent-app
```


- 검색 tool 고도화 예정

<img width="2512" height="1574" alt="image" src="https://github.com/user-attachments/assets/fa1fda74-4278-473e-b5c0-451731fabfe7" />

<img width="2574" height="1399" alt="image" src="https://github.com/user-attachments/assets/307b49cd-1c1f-4fac-9a21-07d217ddac9e" />



