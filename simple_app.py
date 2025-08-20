from fastapi import FastAPI
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Optional
import os
import json
import base64
from io import BytesIO

# HJ ChemAgent 관련 임포트
from langchain_google_genai import ChatGoogleGenerativeAI
from langchain.tools import tool
from typing import Annotated
from langchain_core.messages import BaseMessage, ToolMessage
from typing_extensions import TypedDict
from langgraph.graph import StateGraph, START, END
from langgraph.graph.message import add_messages
import requests
import urllib.parse

app = FastAPI(title="HJ ChemAgent Simple")

# LLM 설정
os.environ["GOOGLE_API_KEY"] = 'Your API KEY'
llm = ChatGoogleGenerativeAI(
    model="gemini-2.0-flash-001",
    temperature=0,
    max_tokens=None,
    timeout=None,
    max_retries=2,
)

# 화학 도구들 정의
@tool
def name_to_smiles_cir(name: str) -> str:
    """IUPAC Name을 SMILES로 변환"""
    encoded_name = urllib.parse.quote(name)
    url = f"https://cactus.nci.nih.gov/chemical/structure/{encoded_name}/smiles"
    
    try:
        response = requests.get(url, timeout=10)
        smiles = response.text.strip()
        return smiles
    except Exception as e:
        return f"Error converting name to SMILES: {str(e)}"

@tool
def generate_molecule_image(smiles: str) -> str:
    """Generate and return base64 encoded image of a molecule from its SMILES string."""
    try:
        # RDKit 임포트 시도
        try:
            from rdkit import Chem
            from rdkit.Chem import Draw
        except ImportError:
            return "RDKit이 설치되지 않았습니다. 분자 이미지를 생성하려면 다음 명령어로 RDKit을 설치해주세요:\n\nconda install -c conda-forge rdkit\n\n또는\n\npip install rdkit"
        
        # SMILES를 분자로 변환
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "Invalid SMILES string."
        
        # 분자 이미지 생성
        image = Draw.MolToImage(mol, size=(300, 300))
        
        # PIL 이미지를 base64로 변환
        buffered = BytesIO()
        image.save(buffered, format="PNG")
        img_base64 = base64.b64encode(buffered.getvalue()).decode()
        
        return f"data:image/png;base64,{img_base64}"
        
    except Exception as e:
        return f"Error generating molecule image: {str(e)}"

tools = [generate_molecule_image, name_to_smiles_cir]

# LangGraph 설정
class State(TypedDict):
    messages: Annotated[list, add_messages]

class BasicToolNode:
    def __init__(self, tools: list) -> None:
        self.tools_by_name = {tool.name: tool for tool in tools}

    def __call__(self, inputs: dict):
        if messages := inputs.get("messages", []):
            message = messages[-1]
        else:
            raise ValueError("No message found in input")
        outputs = []
        for tool_call in message.tool_calls:
            tool_result = self.tools_by_name[tool_call["name"]].invoke(
                tool_call["args"]
            )
            outputs.append(
                ToolMessage(
                    content=json.dumps(tool_result),
                    name=tool_call["name"],
                    tool_call_id=tool_call["id"],
                )
            )
        return {"messages": outputs}

# 그래프 구성
graph_builder = StateGraph(State)
tool_node = BasicToolNode(tools=tools)
graph_builder.add_node("tools", tool_node)

llm_with_tools = llm.bind_tools(tools)

def chatbot(state: State):
    return {"messages": [llm_with_tools.invoke(state["messages"])]}

graph_builder.add_node("chatbot", chatbot)

def route_tools(state: State):
    if isinstance(state, list):
        ai_message = state[-1]
    elif messages := state.get("messages", []):
        ai_message = messages[-1]
    else:
        raise ValueError(f"No messages found in input state to tool_edge: {state}")
    if hasattr(ai_message, "tool_calls") and len(ai_message.tool_calls) > 0:
        return "tools"
    return END

graph_builder.add_conditional_edges(
    "chatbot",
    route_tools,
    {"tools": "tools", END: END},
)
graph_builder.add_edge("tools", "chatbot")
graph_builder.add_edge(START, "chatbot")
graph = graph_builder.compile()

# CORS 설정
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class ChatMessage(BaseModel):
    message: str

class ChatResponse(BaseModel):
    response: str
    image: Optional[str] = None

@app.get("/", response_class=HTMLResponse)
async def get_chat_interface():
    try:
        with open("static/index.html", "r", encoding="utf-8") as f:
            return HTMLResponse(f.read())
    except Exception as e:
        return HTMLResponse(f"<h1>Error loading page: {str(e)}</h1>")

@app.post("/chat")
async def chat_endpoint(chat_message: ChatMessage):
    try:
        # HJ ChemAgent 실행
        result = graph.invoke({"messages": [{"role": "user", "content": chat_message.message}]})
        
        # 마지막 AI 메시지 가져오기
        last_message = result["messages"][-1]
        response_text = last_message.content
        
        # 이미지가 포함된 응답인지 확인
        image_data = None
        for message in result["messages"]:
            if hasattr(message, 'content') and isinstance(message.content, str):
                # JSON으로 인코딩된 툴 결과에서 이미지 데이터 추출
                try:
                    import json
                    content = json.loads(message.content)
                    if isinstance(content, str) and content.startswith("data:image/png;base64,"):
                        image_data = content
                        break
                except:
                    # JSON이 아닌 경우 직접 확인
                    if message.content.startswith("data:image/png;base64,"):
                        image_data = message.content
                        break
        
        # 응답 텍스트에서 base64 이미지 데이터 제거 (텍스트로 표시되지 않도록)
        if image_data and image_data in response_text:
            response_text = response_text.replace(image_data, "").strip()
        
        # base64 데이터가 응답에 포함되어 있다면 제거하고 간단한 메시지로 대체
        if "data:image/png;base64," in response_text:
            response_text = "CCCC 분자의 구조 이미지를 생성했습니다."
        
        return ChatResponse(response=response_text, image=image_data)
        
    except Exception as e:
        return ChatResponse(
            response=f"오류가 발생했습니다: {str(e)}"
        )

# 정적 파일 서빙
try:
    app.mount("/static", StaticFiles(directory="static"), name="static")
except Exception as e:
    print(f"Static files mounting error: {e}")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8002)
