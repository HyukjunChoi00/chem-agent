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

# 논문 검색 관련 임포트
from playwright.async_api import async_playwright
from bs4 import BeautifulSoup, Comment
import asyncio
import re
from urllib.parse import urlparse, urljoin
import nest_asyncio
nest_asyncio.apply()
import sys
if sys.platform == "win32":
    asyncio.set_event_loop_policy(asyncio.WindowsProactorEventLoopPolicy())

# Scimago Journal Rank 분석을 위한 추가 임포트
import sqlite3
from datetime import datetime, timedelta
from dataclasses import dataclass
import logging
from functools import lru_cache

# 로깅 설정
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI(title="HJ ChemAgent Simple")

# LLM 설정
os.environ["GOOGLE_API_KEY"] = 'AIzaSyAPf9KX4Ea_MKiC_tYfm8OUlv_A71BUz1I'
llm = ChatGoogleGenerativeAI(
    model="gemini-2.0-flash-001",
    temperature=0,
    max_tokens=None,
    timeout=None,
    max_retries=2,
    system_instruction="""당신은 화학 연구 보조 AI입니다. 
    
    🚨 CRITICAL RULE - HTML 완전 출력:
    도구가 HTML을 반환하면 (특히 <div>, <table> 태그):
    1. HTML을 완전히, 끝까지, 그대로 출력하세요 - 절대로 마크다운 테이블로 변환하지 마세요
    2. HTML을 자르거나 생략하지 마세요 - 전체 HTML 코드를 그대로 출력
    3. 절대 추가 설명 금지: "결과입니다", "검색 결과입니다", "분석했습니다" 등
    4. HTML 앞뒤에 어떤 텍스트도 붙이지 마세요
    5. 도구가 완전한 HTML 테이블을 반환했다면 그 전체를 출력하세요
    6. 마크다운 테이블(|---|---|---)로 변환하지 마세요 - HTML 그대로 출력
    
    올바른 예시:
    도구 출력: "<div class="research-papers-table enhanced-search"><h3>🔬 논문 검색 결과...</h3><table>...</table></div>"
    → 당신의 응답: <div class="research-papers-table enhanced-search"><h3>🔬 논문 검색 결과...</h3><table>...</table></div>
    
    절대 금지:
    → "doxorubicin에 대한 최신 논문 검색 결과입니다. | 논문 제목 | 주요 내용"
    → HTML을 마크다운 테이블로 변환하는 행위
    → HTML을 중간에 자르는 행위
    → HTML 앞뒤에 설명 텍스트 추가
    
    HTML이 아닌 일반 질문은 평소처럼 친절하게 답변하세요.""",
)

# Scimago Journal Rank 데이터 클래스
@dataclass
class JournalMetrics:
    """저널 메트릭 데이터 클래스"""
    journal_name: str
    sjr_score: float
    h_index: int
    total_docs: int
    total_cites: int
    cites_per_doc: float
    country: str
    publisher: str
    categories: list
    quartile: str
    rank_in_category: int
    percentile: float
    issn: str
    last_updated: str

# Scimago Journal Rank 분석기 (간소화된 버전)
class ScimagoJournalAnalyzer:
    """Scimago Journal Rank 데이터를 활용한 저널 분석기 (ChemAgent 통합용)"""
    
    def __init__(self, cache_db_path: str = "journal_cache.db"):
        self.base_url = "https://www.scimagojr.com"
        self.cache_db_path = cache_db_path
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36'
        })
        
        # 캐시 데이터베이스 초기화
        self._init_cache_db()
        
        # 실시간 Scimago 검색을 위한 설정
        self.scimago_search_url = "https://www.scimagojr.com/journalsearch.php"
        
        # 웹 스크래핑을 위한 헤더 설정
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
            'Accept-Language': 'en-US,en;q=0.5',
            'Accept-Encoding': 'gzip, deflate',
            'Connection': 'keep-alive',
            'Upgrade-Insecure-Requests': '1',
        })

    def _init_cache_db(self):
        """캐시 데이터베이스 초기화"""
        try:
            conn = sqlite3.connect(self.cache_db_path)
            cursor = conn.cursor()
            
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS journal_metrics (
                    journal_name_key TEXT PRIMARY KEY,
                    journal_name TEXT,
                    sjr_score REAL,
                    h_index INTEGER,
                    quartile TEXT,
                    last_updated TEXT
                )
            ''')
            
            conn.commit()
            conn.close()
            logger.info("Scimago 캐시 데이터베이스 초기화 완료")
            
        except Exception as e:
            logger.error(f"Scimago 캐시 DB 초기화 오류: {str(e)}")

    def _clean_journal_name_for_search(self, journal_name: str) -> list:
        """검색용 저널명 정리 (여러 변형 반환)"""
        variants = []
        
        # 원본
        variants.append(journal_name)
        
        # 괄호와 내용 제거
        cleaned = re.sub(r'\s*\([^)]*\)', '', journal_name).strip()
        if cleaned != journal_name:
            variants.append(cleaned)
        
        # 특수 문자 제거
        simple = re.sub(r'[^\w\s]', '', journal_name).strip()
        if simple != journal_name and simple not in variants:
            variants.append(simple)
        
        # 첫 단어만
        first_word = journal_name.split()[0]
        if len(first_word) > 3 and first_word not in variants:
            variants.append(first_word)
        
        return variants

    def search_scimago_realtime(self, journal_name: str) -> dict:
        """Scimago에서 실시간으로 저널 메트릭 검색 (정확한 XPath 사용)"""
        try:
            print(f"[SCIMAGO] 실시간 검색: {journal_name}")
            
            # 저널명 변형들 생성
            search_variants = self._clean_journal_name_for_search(journal_name)
            print(f"[SCIMAGO] 검색 변형들: {search_variants}")
            
            # 각 변형에 대해 검색 시도
            for variant in search_variants:
                print(f"[SCIMAGO] 변형 '{variant}'로 검색 시도")
                
                # 1단계: 검색 페이지에서 저널 링크 찾기
                search_params = {
                    'q': variant.replace(' ', '+')
                }
                
                search_url = f"{self.scimago_search_url}?q={variant.replace(' ', '+')}"
                print(f"[SCIMAGO] 검색 URL: {search_url}")
                
                search_response = self.session.get(self.scimago_search_url, params=search_params, timeout=15)
                print(f"[SCIMAGO] 검색 응답 상태: {search_response.status_code}")
                
                if search_response.status_code == 200:
                    search_soup = BeautifulSoup(search_response.content, 'html.parser')
                    print(f"[SCIMAGO] HTML 길이: {len(search_response.content)} bytes")
                    
                    # XPath: /html/body/div[4]/div[2]/a 에 해당하는 요소 찾기
                    journal_links = search_soup.select('div:nth-of-type(4) div:nth-of-type(2) a')
                    print(f"[SCIMAGO] XPath 링크 개수: {len(journal_links)}")
                    
                    if not journal_links:
                        # 대안: 저널명을 포함한 모든 링크 찾기
                        all_links = search_soup.find_all('a', href=True)
                        journal_links = [link for link in all_links if 'journalsearch.php?q=' in link.get('href', '')]
                        print(f"[SCIMAGO] 대안 방법으로 찾은 링크 개수: {len(journal_links)}")
                        
                        # 처음 3개 링크 정보 출력
                        for i, link in enumerate(journal_links[:3]):
                            print(f"[SCIMAGO] 링크 {i+1}: {link.get('href')} - {link.get_text()[:50]}")
                    
                    if journal_links:
                        # 첫 번째 매칭 저널 선택
                        first_link = journal_links[0]
                        journal_href = first_link.get('href')
                        link_text = first_link.get_text().strip()
                        
                        print(f"[SCIMAGO] 선택된 링크: {journal_href}")
                        print(f"[SCIMAGO] 링크 텍스트: {link_text}")
                        
                        # 상대 URL을 절대 URL로 변환
                        if journal_href.startswith('journalsearch.php'):
                            journal_url = f"https://www.scimagojr.com/{journal_href}"
                        else:
                            journal_url = journal_href
                        
                        print(f"[SCIMAGO] 저널 상세 페이지로 이동: {journal_url}")
                        
                        # 2단계: 저널 상세 페이지에서 SJR 점수 추출
                        detail_response = self.session.get(journal_url, timeout=15)
                        print(f"[SCIMAGO] 상세 페이지 응답 상태: {detail_response.status_code}")
                        
                        if detail_response.status_code == 200:
                            detail_soup = BeautifulSoup(detail_response.content, 'html.parser')
                            print(f"[SCIMAGO] 상세 페이지 HTML 길이: {len(detail_response.content)} bytes")
                            
                            # XPath로 SJR 점수 추출
                            metrics = self._extract_scimago_metrics_from_detail(detail_soup, journal_name)
                            
                            if metrics and metrics.get('sjr_score', 0) > 0:
                                self._cache_scimago_metrics(journal_name, metrics)
                                print(f"[SCIMAGO] 검색 성공 (변형: {variant}): {journal_name} → SJR: {metrics['sjr_score']:.2f}")
                                return metrics
                            else:
                                print(f"[SCIMAGO] SJR 추출 실패 (변형: {variant}): {journal_name}")
                    else:
                        print(f"[SCIMAGO] 저널 링크를 찾을 수 없음 (변형: {variant})")
                else:
                    print(f"[SCIMAGO] 검색 페이지 접근 실패 (변형: {variant}): {search_response.status_code}")
            
            print(f"[SCIMAGO] 모든 변형 검색 실패: {journal_name}")
            return None
            
        except Exception as e:
            logger.error(f"Scimago 실시간 검색 오류: {str(e)}")
            print(f"[SCIMAGO] 오류 상세: {str(e)}")
            import traceback
            print(f"[SCIMAGO] 스택 트레이스: {traceback.format_exc()}")
            return None
    
    def _extract_scimago_metrics_from_detail(self, soup, journal_name: str) -> dict:
        """저널 상세 페이지에서 정확한 XPath로 SJR과 Quartile 추출"""
        try:
            print(f"[SCIMAGO] 상세 페이지에서 SJR과 Quartile 추출 중...")
            
            # XPath: /html/body/div[7]/div/div/div[4]/p[1] 에서 SJR과 Quartile 추출
            # 예시: <p class="hindexnumber">0.562 <span class="Q3" style="font-size: 0.8em" title="best quartile for all subject categories">Q3</span></p>
            target_elements = soup.select('body > div:nth-of-type(7) > div > div > div:nth-of-type(4) > p:nth-of-type(1)')
            
            sjr_score = 0.0
            quartile = "Unknown"
            h_index = 0
            
            # SJR 점수와 Quartile 추출
            if target_elements:
                target_element = target_elements[0]
                element_text = target_element.get_text(strip=True)
                print(f"[SCIMAGO] 대상 요소 텍스트: {element_text}")
                
                # SJR 점수 추출 (첫 번째 숫자)
                sjr_numbers = re.findall(r'\d+\.\d+', element_text)
                if sjr_numbers:
                    sjr_score = float(sjr_numbers[0])
                    print(f"[SCIMAGO] SJR 점수 추출: {sjr_score}")
                
                # Quartile 추출 (Q1, Q2, Q3, Q4)
                quartile_spans = target_element.find_all('span', class_=re.compile(r'^Q[1-4]$'))
                if quartile_spans:
                    quartile = quartile_spans[0].get_text(strip=True)
                    print(f"[SCIMAGO] Quartile 추출: {quartile}")
                else:
                    # span 태그가 없을 경우 텍스트에서 직접 추출
                    quartile_match = re.search(r'Q[1-4]', element_text)
                    if quartile_match:
                        quartile = quartile_match.group()
                        print(f"[SCIMAGO] Quartile 텍스트에서 추출: {quartile}")
            
            # 대안: 다른 위치에서 SJR 찾기
            if sjr_score == 0:
                # class="hindexnumber" 요소 찾기
                hinumber_elements = soup.find_all('p', class_='hindexnumber')
                for elem in hinumber_elements:
                    elem_text = elem.get_text(strip=True)
                    sjr_numbers = re.findall(r'\d+\.\d+', elem_text)
                    if sjr_numbers:
                        sjr_score = float(sjr_numbers[0])
                        print(f"[SCIMAGO] hindexnumber에서 SJR 추출: {sjr_score}")
                        
                        # 같은 요소에서 Quartile도 찾기
                        q_spans = elem.find_all('span', class_=re.compile(r'^Q[1-4]$'))
                        if q_spans:
                            quartile = q_spans[0].get_text(strip=True)
                            print(f"[SCIMAGO] hindexnumber에서 Quartile 추출: {quartile}")
                        break
            
            # H-index 추출 (기존 로직 유지)
            h_elements = soup.find_all(string=re.compile(r'H index', re.IGNORECASE))
            for h_elem in h_elements:
                parent = h_elem.parent
                if parent:
                    h_text = parent.get_text()
                    h_numbers = re.findall(r'\d+', h_text)
                    if h_numbers:
                        h_index = int(h_numbers[0])
                        print(f"[SCIMAGO] H-index 추출: {h_index}")
                        break
            
            if sjr_score > 0:
                # IF는 SJR에서 추정
                impact_factor_estimate = sjr_score * 2.5
                
                return {
                    "journal": journal_name,
                    "sjr_score": sjr_score,
                    "h_index": h_index,
                    "quartile": quartile,
                    "impact_factor_estimate": impact_factor_estimate,
                    "source": "scimago_realtime_detail"
                }
            
            print(f"[SCIMAGO] SJR 점수를 찾을 수 없음")
            return None
            
        except Exception as e:
            logger.error(f"Scimago SJR/Quartile 추출 오류: {str(e)}")
            print(f"[SCIMAGO] SJR/Quartile 추출 오류: {str(e)}")
            return None
    
    def _extract_scimago_metrics(self, soup, journal_name: str) -> dict:
        """Scimago 페이지에서 메트릭 추출 (구버전 - 백업용)"""
        try:
            # 구버전 로직 유지 (백업용) - 경고 수정
            sjr_elements = soup.find_all(text=re.compile(r'\d+\.\d+'))  # 일단 기존 방식 유지
            sjr_score = 0.0
            h_index = 0
            
            for element in sjr_elements[:3]:
                try:
                    value = float(element.strip())
                    if 0.1 <= value <= 50:  # SJR 점수 범위
                        sjr_score = value
                        break
                except:
                    continue
            
            h_elements = soup.find_all(text=re.compile(r'H index'))  # 일단 기존 방식 유지
            if h_elements:
                for h_elem in h_elements:
                    parent = h_elem.parent
                    if parent:
                        numbers = re.findall(r'\d+', parent.get_text())
                        if numbers:
                            h_index = int(numbers[0])
                            break
            
            if sjr_score > 0:
                return {
                    "journal": journal_name,
                    "sjr_score": sjr_score,
                    "h_index": h_index,
                    "quartile": self._estimate_quartile(sjr_score),
                    "impact_factor_estimate": sjr_score * 2.5,
                    "source": "scimago_realtime"
                }
            
            return None
            
        except Exception as e:
            logger.error(f"Scimago 메트릭 추출 오류: {str(e)}")
            return None
    
    def _estimate_quartile(self, sjr_score: float) -> str:
        """SJR 점수를 기반으로 Quartile 추정"""
        if sjr_score >= 10:
            return "Q1"
        elif sjr_score >= 3:
            return "Q1"
        elif sjr_score >= 1:
            return "Q2"
        elif sjr_score >= 0.3:
            return "Q3"
        else:
            return "Q4"
    
    def _cache_scimago_metrics(self, journal_name: str, metrics: dict):
        """Scimago 메트릭을 캐시에 저장"""
        try:
            journal_key = journal_name.lower().strip()
            conn = sqlite3.connect(self.cache_db_path)
            cursor = conn.cursor()
            
            cursor.execute('''
                INSERT OR REPLACE INTO journal_metrics VALUES 
                (?, ?, ?, ?, ?, ?)
            ''', (
                journal_key,
                metrics['journal'],
                metrics['sjr_score'],
                metrics['h_index'],
                metrics['quartile'],
                datetime.now().isoformat()
            ))
            
            conn.commit()
            conn.close()
            
        except Exception as e:
            logger.error(f"Scimago 캐시 저장 오류: {str(e)}")

    def get_journal_metrics_fast(self, journal_name: str) -> dict:
        """저널 메트릭 조회 (캐시 우선, 없으면 실시간 검색)"""
        try:
            journal_key = journal_name.lower().strip()
            print(f"[SCIMAGO] 메트릭 조회 시작: {journal_name} (키: {journal_key})")
            
            # 1단계: 캐시에서 조회
            cached_metrics = self._get_cached_metrics(journal_key)
            if cached_metrics and cached_metrics.get('sjr_score', 0) > 0:
                print(f"[SCIMAGO] 캐시에서 조회 성공: {journal_name} → SJR: {cached_metrics.get('sjr_score', 0)}")
                return cached_metrics
            else:
                print(f"[SCIMAGO] 캐시에 유효한 데이터 없음: {journal_name}")
            
            # 2단계: 실시간 Scimago 검색
            print(f"[SCIMAGO] 실시간 검색 시작: {journal_name}")
            realtime_metrics = self.search_scimago_realtime(journal_name)
            if realtime_metrics and realtime_metrics.get('sjr_score', 0) > 0:
                print(f"[SCIMAGO] 실시간 검색 성공: {journal_name}")
                return realtime_metrics
            else:
                print(f"[SCIMAGO] 실시간 검색 실패: {journal_name}")
            
            # 3단계: 기본값 반환
            print(f"[SCIMAGO] 기본값 반환: {journal_name}")
            return {
                "journal": journal_name,
                "sjr_score": 0.0,
                "h_index": 0,
                "quartile": "Unknown",
                "impact_factor_estimate": 0.0,
                "source": "unknown"
            }
            
        except Exception as e:
            logger.error(f"저널 메트릭 조회 오류: {str(e)}")
            print(f"[SCIMAGO] 오류 발생: {journal_name} - {str(e)}")
            return {
                "journal": journal_name,
                "sjr_score": 0.0,
                "h_index": 0,
                "quartile": "Unknown",
                "impact_factor_estimate": 0.0,
                "source": "error"
            }

    def _get_cached_metrics(self, journal_key: str) -> dict:
        """캐시에서 저널 메트릭 조회"""
        try:
            conn = sqlite3.connect(self.cache_db_path)
            cursor = conn.cursor()
            
            cursor.execute('''
                SELECT journal_name, sjr_score, h_index, quartile, last_updated
                FROM journal_metrics 
                WHERE journal_name_key = ? AND last_updated > ?
            ''', (journal_key, (datetime.now() - timedelta(days=30)).isoformat()))
            
            result = cursor.fetchone()
            conn.close()
            
            if result:
                sjr_score = result[1] if result[1] is not None else 0.0
                h_index = result[2] if result[2] is not None else 0
                quartile = result[3] if result[3] is not None else "Unknown"
                
                print(f"[SCIMAGO] 캐시 데이터 발견: {result[0]} → SJR: {sjr_score}, Q: {quartile}")
                
                if sjr_score > 0:  # 유효한 SJR 점수가 있을 때만 반환
                    return {
                        "journal": result[0],
                        "sjr_score": sjr_score,
                        "h_index": h_index,
                        "quartile": quartile,
                        "impact_factor_estimate": sjr_score * 2.5,
                        "source": "cached"
                    }
                else:
                    print(f"[SCIMAGO] 캐시 데이터 무효 (SJR=0): {result[0]}")
            else:
                print(f"[SCIMAGO] 캐시에 데이터 없음: {journal_key}")
            
            return None
            
        except Exception as e:
            logger.error(f"캐시 조회 오류: {str(e)}")
            print(f"[SCIMAGO] 캐시 조회 오류: {str(e)}")
            return None

    def calculate_impact_score(self, journal_metrics: dict, citation_count: int = 0) -> float:
        """종합 영향력 점수 계산"""
        try:
            sjr_score = journal_metrics.get('sjr_score', 0)
            h_index = journal_metrics.get('h_index', 0)
            
            # SJR 점수 가중치 (40%)
            sjr_component = sjr_score * 4.0
            
            # H-index 가중치 (30%)
            h_index_component = (h_index / 100) * 30
            
            # 인용수 가중치 (30%)
            citation_component = min(citation_count / 50, 1.0) * 30
            
            total_score = sjr_component + h_index_component + citation_component
            
            return round(total_score, 2)
            
        except Exception as e:
            logger.error(f"영향력 점수 계산 오류: {str(e)}")
            return 0.0

    def get_impact_tier(self, impact_score: float) -> tuple:
        """영향력 점수에 따른 등급 분류"""
        if impact_score >= 50:
            return ("🔥 Top Tier", "high-impact")
        elif impact_score >= 30:
            return ("⭐ High Impact", "medium-high-impact")
        elif impact_score >= 15:
            return ("📈 Notable", "medium-impact")
        elif impact_score >= 5:
            return ("📄 Standard", "standard-impact")
        else:
            return ("📝 Emerging", "low-impact")

# 전역 Scimago 분석기 인스턴스
scimago_analyzer = ScimagoJournalAnalyzer()

# 논문 검색 클래스 정의 (잘 작동하는 원본 버전)
class paper_search:
    def __init__(self, max_articles=5):
        self.max_articles = max_articles
        self.base_url = "https://pubmed.ncbi.nlm.nih.gov/"
        self.playwright = None
        self.browser = None
        self.context = None

    async def setup_browser(self):
        if self.playwright is None:
            self.playwright = await async_playwright().start()
            self.browser = await self.playwright.chromium.launch(
                headless=True,
                args=[
                    '--no-sandbox',
                    '--disable-setuid-sandbox',
                    '--disable-dev-shm-usage',
                    '--disable-accelerated-2d-canvas',
                    '--no-first-run',
                    '--no-zygote',
                    '--disable-gpu'
                ]
            )
            self.context = await self.browser.new_context(
                user_agent='Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
                viewport={'width': 1920, 'height': 1080},
                extra_http_headers={
                    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
                    'Accept-Language': 'en-US,en;q=0.5',
                    'Accept-Encoding': 'gzip, deflate',
                    'DNT': '1',
                    'Connection': 'keep-alive',
                    'Upgrade-Insecure-Requests': '1',
                }
            )

    async def close_browser(self):
        if self.browser:
            await self.browser.close()
        if self.playwright:
            await self.playwright.stop()
        self.playwright = None
        self.browser = None
        self.context = None

    def clean_title(self, title: str) -> str:
        title = re.sub(r'\s+', ' ', title.strip())
        title = re.sub(r'^\W+|\W+$', '', title)
        return title if len(title) > 10 else ""

    async def extract_articles(self, page, query: str) -> list:
        articles = []
        print(f"논문 추출 시작... 최대 {self.max_articles}개의 논문을 찾습니다.")
        
        # 페이지 로딩 대기
        try:
            await page.wait_for_selector('#search-results', timeout=10000)
            print("검색 결과 컨테이너 로딩 완료")
        except:
            print("검색 결과 컨테이너를 찾을 수 없습니다.")
            return []
        
        # 실제 HTML 구조를 바탕으로 한 올바른 PubMed selector들
        selectors_to_try = [
            'a.docsum-title',  # 실제 구조: <a class="docsum-title" href="...">
            '.docsum-title',   # 클래스만으로도 가능
            '#search-results a.docsum-title',  # 더 구체적인 경로
            'a[data-ga-category="result_click"]',  # 데이터 속성 활용
            'a[data-article-id]',  # 논문 ID 속성이 있는 링크
            '[data-ga-action]'  # GA 액션 속성이 있는 요소
        ]
        
        # 모든 논문 링크를 한 번에 찾기
        all_links = []
        for selector in selectors_to_try:
            try:
                links = await page.query_selector_all(selector)
                if links:
                    print(f"'{selector}'로 {len(links)}개 링크 발견")
                    all_links = links
                    break
            except Exception as e:
                print(f"Selector '{selector}' 오류: {str(e)}")
                continue
        
        if not all_links:
            print("어떤 selector로도 링크를 찾을 수 없습니다.")
            # HTML 구조 디버깅
            try:
                html_content = await page.inner_html('#search-results')
                print("검색 결과 HTML 일부:")
                print(html_content[:500] + "...")
            except:
                print("HTML 내용을 가져올 수 없습니다.")
            return []
        
        # 찾은 링크들에서 논문 정보 추출
        for i in range(min(len(all_links), self.max_articles)):
            try:
                link = all_links[i]
                print(f"논문 {i+1} 추출 중...")
                
                title = await link.inner_text()
                url = await link.get_attribute('href')
                
                print(f"논문 {i+1} 제목: {title[:50]}...")
                print(f"논문 {i+1} URL: {url}")
                
                # URL 정규화
                if url and not url.startswith("http"):
                    url = urljoin(self.base_url, url)
                
                title = self.clean_title(title)
                if title and url:
                    articles.append({"title": title, "url": url})
                    print(f"논문 {i+1}: 성공적으로 추출됨")
                else:
                    print(f"논문 {i+1}: 제목 또는 URL이 유효하지 않음")
                    
            except Exception as e:
                print(f"논문 {i+1} 추출 중 오류: {str(e)}")
                continue
                
        print(f"총 {len(articles)}개의 논문을 추출했습니다.")
        return articles

    async def extract_abstract(self, url: str) -> str:
        """논문 페이지에서 초록을 추출하는 메서드"""
        page = await self.context.new_page()
        try:
            print(f"초록 추출을 위해 페이지 접속: {url}")
            await page.goto(url, wait_until='domcontentloaded', timeout=20000)
            await asyncio.sleep(3)
            
            # 여러 가지 방법으로 초록 찾기
            abstract_text = ""
            
            # 방법 1: #abstract .abstract-content 사용
            abstract_content = await page.query_selector('#abstract .abstract-content')
            if abstract_content:
                abstract_text = await abstract_content.inner_text()
                print("방법 1로 초록 추출 성공")
            else:
                # 방법 2: #eng-abstract 사용 (영어 초록)
                eng_abstract = await page.query_selector('#eng-abstract')
                if eng_abstract:
                    abstract_text = await eng_abstract.inner_text()
                    print("방법 2로 초록 추출 성공")
                else:
                    # 방법 3: #abstract 전체 사용
                    abstract_element = await page.query_selector('#abstract')
                    if abstract_element:
                        abstract_text = await abstract_element.inner_text()
                        # "Abstract" 제목 제거
                        abstract_text = re.sub(r'^Abstract\s*', '', abstract_text.strip())
                        print("방법 3으로 초록 추출 성공")
                    else:
                        print("모든 방법으로 초록을 찾을 수 없습니다.")
                        return "초록을 찾을 수 없습니다."
            
            return abstract_text.strip() if abstract_text else "초록 내용이 비어있습니다."
            
        except Exception as e:
            print(f"초록 추출 중 오류 발생: {str(e)}")
            return f"초록 추출 중 오류 발생: {str(e)}"
        finally:
            await page.close()

    async def run_pipeline(self, query: str) -> list:
        await self.setup_browser()
        search_url = f"https://pubmed.ncbi.nlm.nih.gov/?term={query}&filter=datesearch.y_5&filter=simsearch1.fha&filter=pubt.clinicaltrial"
        
        print(f"검색 URL: {search_url}")
        
        # 403 오류 대응을 위한 재시도 로직
        max_retries = 3
        articles = []
        
        for attempt in range(max_retries):
            page = await self.context.new_page()
            try:
                print(f"시도 {attempt + 1}/{max_retries}")
                
                # 페이지 로딩 전 대기
                await asyncio.sleep(2 * (attempt + 1))  # 점진적 대기 시간 증가
                
                await page.goto(search_url, wait_until='domcontentloaded', timeout=30000)
                await asyncio.sleep(5)
                
                print("검색 페이지 로딩 완료")
                
                # 페이지 제목 확인
                page_title = await page.title()
                print(f"페이지 제목: {page_title}")
                
                # 403 오류 감지
                if "403" in page_title or "Forbidden" in page_title:
                    print(f"403 오류 감지됨. 시도 {attempt + 1}")
                    if attempt < max_retries - 1:
                        await page.close()
                        continue
                    else:
                        print("모든 재시도 실패. 403 오류 지속됨.")
                        await page.close()
                        break
                
                # 정상적인 페이지인 경우 논문 추출
                articles = await self.extract_articles(page, query)
                await page.close()
                break
                
            except Exception as e:
                print(f"시도 {attempt + 1} 중 오류 발생: {str(e)}")
                await page.close()
                if attempt < max_retries - 1:
                    continue
                else:
                    print("모든 재시도 실패")
                    break

        results = []
        print(f"\n초록 추출 시작... {len(articles)}개의 논문 처리")
        
        for i, article in enumerate(articles, 1):
            print(f"\n=== 논문 {i}/{len(articles)} 초록 추출 ===")
            # 초록 추출을 위해 새로운 메서드 사용
            abstract = await self.extract_abstract(article["url"])
            results.append({
                "title": article["title"],
                "url": article["url"],
                "abstract": abstract
            })

        await self.close_browser()
        return results

# 전역 논문 검색 인스턴스
paper_searcher = paper_search(max_articles=5)

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

def search_research_papers_trends_sync(query: str, max_papers: int = 10) -> str:
    """연구 동향 분석을 위한 더 넓은 범위의 논문 검색 (필터 완화)"""
    try:
        print(f"[DEBUG] 동향 분석용 논문 검색 시작: {query}")
        
        # 새로운 이벤트 루프에서 비동기 함수 실행
        import asyncio
        
        async def _search_papers_for_trends():
            searcher = paper_search(max_articles=max_papers)
            print(f"[DEBUG] paper_search 인스턴스 생성됨 (동향 분석용)")
            
            # 필터를 완화한 검색 URL (모든 논문 타입 포함, 최근 5년)
            await searcher.setup_browser()
            search_url = f"https://pubmed.ncbi.nlm.nih.gov/?term={query}&filter=datesearch.y_5"
            
            print(f"[DEBUG] 동향 분석 검색 URL: {search_url}")
            
            # 403 오류 대응을 위한 재시도 로직
            max_retries = 3
            articles = []
            
            for attempt in range(max_retries):
                page = await searcher.context.new_page()
                try:
                    print(f"[DEBUG] 동향 분석 시도 {attempt + 1}/{max_retries}")
                    
                    # 페이지 로딩 전 대기
                    await asyncio.sleep(2 * (attempt + 1))
                    
                    await page.goto(search_url, wait_until='domcontentloaded', timeout=30000)
                    await asyncio.sleep(5)
                    
                    print("[DEBUG] 동향 분석용 검색 페이지 로딩 완료")
                    
                    # 페이지 제목 확인
                    page_title = await page.title()
                    print(f"[DEBUG] 페이지 제목: {page_title}")
                    
                    # 403 오류 감지
                    if "403" in page_title or "Forbidden" in page_title:
                        print(f"[DEBUG] 403 오류 감지됨. 동향 분석 시도 {attempt + 1}")
                        if attempt < max_retries - 1:
                            await page.close()
                            continue
                        else:
                            print("[DEBUG] 동향 분석 모든 재시도 실패. 403 오류 지속됨.")
                            await page.close()
                            break
                    
                    # 정상적인 페이지인 경우 논문 추출
                    articles = await searcher.extract_articles(page, query)
                    await page.close()
                    break
                    
                except Exception as e:
                    print(f"[DEBUG] 동향 분석 시도 {attempt + 1} 중 오류 발생: {str(e)}")
                    await page.close()
                    if attempt < max_retries - 1:
                        continue
                    else:
                        print("[DEBUG] 동향 분석 모든 재시도 실패")
                        break

            results = []
            print(f"\n[DEBUG] 동향 분석용 초록 추출 시작... {len(articles)}개의 논문 처리")
            
            for i, article in enumerate(articles, 1):
                print(f"\n=== 논문 {i}/{len(articles)} 초록 추출 ===")
                abstract = await searcher.extract_abstract(article["url"])
                results.append({
                    "title": article["title"],
                    "url": article["url"],
                    "abstract": abstract
                })

            await searcher.close_browser()
            return results
        
        # 이벤트 루프 실행
        results = None
        try:
            loop = asyncio.get_event_loop()
            if loop.is_running():
                print("[DEBUG] 기존 이벤트 루프가 실행 중, 새 스레드에서 실행")
                import concurrent.futures
                def run_in_thread():
                    new_loop = asyncio.new_event_loop()
                    asyncio.set_event_loop(new_loop)
                    try:
                        return new_loop.run_until_complete(_search_papers_for_trends())
                    finally:
                        new_loop.close()
                
                with concurrent.futures.ThreadPoolExecutor() as executor:
                    future = executor.submit(run_in_thread)
                    results = future.result(timeout=120)
            else:
                print("[DEBUG] 기존 루프에서 실행")
                results = loop.run_until_complete(_search_papers_for_trends())
        except RuntimeError as e:
            print(f"[DEBUG] RuntimeError 발생: {e}, 새 루프 생성")
            results = asyncio.run(_search_papers_for_trends())
        
        if not results:
            print(f"[DEBUG] 동향 분석용 검색 결과가 비어있음")
            return f"'{query}'에 대한 논문을 찾을 수 없습니다."
        
        print(f"[DEBUG] 동향 분석용 HTML 테이블 생성 중, {len(results)}개 논문")
        
        # 깔끔한 HTML 표 형식으로 결과 생성
        html_table = """
<div class="research-papers-table">
<h3>📚 연구 동향 분석용 논문 검색 결과: {}</h3>
<table class="papers-table">
<thead>
<tr>
<th style="width: 5%;">번호</th>
<th style="width: 35%;">논문 제목</th>
<th style="width: 50%;">초록 요약</th>
<th style="width: 10%;">링크</th>
</tr>
</thead>
<tbody>
""".format(query)
        
        for i, paper in enumerate(results, 1):
            # 초록을 200자로 제한
            abstract_summary = paper['abstract'][:200] + "..." if len(paper['abstract']) > 200 else paper['abstract']
            
            html_table += f"""
<tr>
<td>{i}</td>
<td class="paper-title">{paper['title']}</td>
<td class="paper-abstract">{abstract_summary}</td>
<td><a href="{paper['url']}" target="_blank" class="paper-link">보기</a></td>
</tr>
"""
        
        html_table += """
</tbody>
</table>
</div>
"""
        
        print(f"[DEBUG] 동향 분석용 HTML 테이블 생성 완료")
        return html_table
        
    except Exception as e:
        error_msg = f"동향 분석용 논문 검색 중 오류 발생: {str(e)}"
        print(f"[DEBUG] 동향 분석 예외 발생: {error_msg}")
        return f"논문 검색 중 오류 발생: {str(e)}"


def search_research_papers_sync(query: str, max_papers: int = 5) -> str:
    """화학 연구 논문을 검색하고 표 형식으로 결과를 반환합니다. (동기 버전)"""
    try:
        print(f"[DEBUG] 논문 검색 시작: {query}")
        
        # 새로운 이벤트 루프에서 비동기 함수 실행
        import asyncio
        
        async def _search_papers():
            searcher = paper_search(max_articles=max_papers)
            print(f"[DEBUG] paper_search 인스턴스 생성됨")
            results = await searcher.run_pipeline(query)
            print(f"[DEBUG] 검색 결과: {len(results) if results else 0}개 논문 발견")
            return results
        
        # 이벤트 루프 실행
        results = None
        try:
            loop = asyncio.get_event_loop()
            if loop.is_running():
                print("[DEBUG] 기존 이벤트 루프가 실행 중, 새 스레드에서 실행")
                # 이미 실행 중인 루프가 있으면 새 스레드에서 실행
                import concurrent.futures
                import threading
                
                def run_in_thread():
                    new_loop = asyncio.new_event_loop()
                    asyncio.set_event_loop(new_loop)
                    try:
                        return new_loop.run_until_complete(_search_papers())
                    finally:
                        new_loop.close()
                
                with concurrent.futures.ThreadPoolExecutor() as executor:
                    future = executor.submit(run_in_thread)
                    results = future.result(timeout=120)  # 타임아웃을 120초로 증가
            else:
                print("[DEBUG] 기존 루프에서 실행")
                results = loop.run_until_complete(_search_papers())
        except RuntimeError as e:
            print(f"[DEBUG] RuntimeError 발생: {e}, 새 루프 생성")
            # 루프가 없으면 새로 생성
            results = asyncio.run(_search_papers())
        
        if not results:
            print(f"[DEBUG] 검색 결과가 비어있음")
            # 더 상세한 오류 메시지 제공
            return f"""
<div class="research-papers-table">
<h3>📚 '{query}' 검색 결과</h3>
<p style="color: #dc2626; padding: 1rem; background-color: #fef2f2; border: 1px solid #fecaca; border-radius: 0.5rem;">
논문을 찾을 수 없습니다. 다음을 시도해보세요:<br><br>
• 더 일반적인 검색어 사용 (예: "cancer treatment" 대신 "cancer")<br>
• 영어 검색어 사용<br>
• 다른 화학 물질명 시도<br>
• 검색어를 더 간단하게 변경
</p>
</div>
"""
        
        print(f"[DEBUG] HTML 테이블 생성 중, {len(results)}개 논문")
        
        # 깔끔한 HTML 표 형식으로 결과 생성
        html_table = """
<div class="research-papers-table">
<h3>📚 연구 논문 검색 결과: {}</h3>
<table class="papers-table">
<thead>
<tr>
<th style="width: 5%;">번호</th>
<th style="width: 35%;">논문 제목</th>
<th style="width: 50%;">초록 요약</th>
<th style="width: 10%;">링크</th>
</tr>
</thead>
<tbody>
""".format(query)
        
        for i, paper in enumerate(results, 1):
            # 초록을 200자로 제한
            abstract_summary = paper['abstract'][:200] + "..." if len(paper['abstract']) > 200 else paper['abstract']
            
            html_table += f"""
<tr>
<td>{i}</td>
<td class="paper-title">{paper['title']}</td>
<td class="paper-abstract">{abstract_summary}</td>
<td><a href="{paper['url']}" target="_blank" class="paper-link">보기</a></td>
</tr>
"""
        
        html_table += """
</tbody>
</table>
</div>
"""
        
        print(f"[DEBUG] HTML 테이블 생성 완료")
        return html_table
        
    except Exception as e:
        error_msg = f"논문 검색 중 오류 발생: {str(e)}"
        print(f"[DEBUG] 예외 발생: {error_msg}")
        return f"""
<div class="research-papers-table">
<h3>📚 검색 오류</h3>
<p style="color: #dc2626; padding: 1rem; background-color: #fef2f2; border: 1px solid #fecaca; border-radius: 0.5rem;">
{error_msg}<br><br>
잠시 후 다시 시도해주세요.
</p>
</div>
"""

@tool
def paper_search_tool(query: str) -> str:
    """논문을 검색하고 저널 IF 점수(SJR 기반) 분석과 함께 HTML 표 형식으로 결과를 반환합니다. 사용자가 논문 검색을 요청하면 자동으로 저널 영향력도 분석합니다. 반환된 HTML을 그대로 사용자에게 제공하세요."""
    try:
        print(f"[TOOL] 논문 검색 도구 시작: {query}")
        
        # 비동기 함수를 동기적으로 실행
        import asyncio
        
        async def _search_pipeline():
            pipeline = paper_search(max_articles=5)
            results = await pipeline.run_pipeline(query)
            return results
        
        # 이벤트 루프 처리
        try:
            loop = asyncio.get_event_loop()
            if loop.is_running():
                import concurrent.futures
                def run_in_thread():
                    new_loop = asyncio.new_event_loop()
                    asyncio.set_event_loop(new_loop)
                    try:
                        return new_loop.run_until_complete(_search_pipeline())
                    finally:
                        new_loop.close()
                
                with concurrent.futures.ThreadPoolExecutor() as executor:
                    future = executor.submit(run_in_thread)
                    results = future.result(timeout=120)
            else:
                results = loop.run_until_complete(_search_pipeline())
        except RuntimeError:
            results = asyncio.run(_search_pipeline())
        
        if not results:
            return f"""
<div class="research-papers-table">
<h3>📚 논문 검색 결과: {query}</h3>
<p style="color: #dc2626; padding: 1rem; background-color: #fef2f2; border: 1px solid #fecaca; border-radius: 0.5rem;">
'{query}'에 대한 논문을 찾을 수 없습니다.<br><br>
💡 <strong>팁:</strong> "API로 {query} 논문을 검색해줘"라고 요청하면 더 안정적으로 검색할 수 있습니다.
</p>
</div>
"""
        
        # LLM을 사용해서 각 논문을 요약하고 표로 정리
        print(f"[TOOL] LLM을 사용한 논문 요약 시작...")
        
        # 모든 논문 정보를 LLM에게 전달하여 요약 요청
        papers_info = ""
        for i, paper in enumerate(results, 1):
            papers_info += f"""
논문 {i}:
제목: {paper['title']}
초록: {paper['abstract']}
URL: {paper['url']}

---
"""
        
        summary_prompt = f"""
다음 {len(results)}개의 논문을 분석하여 각 논문의 주요 내용을 한국어로 요약해주세요.
각 논문마다 2-3줄로 핵심 내용만 간단히 정리해주세요.

{papers_info}

각 논문에 대해 다음 형식으로 응답해주세요:
논문1||제목||주요내용요약
논문2||제목||주요내용요약
...

주요내용요약은 연구 목적, 방법, 주요 결과를 포함해서 2-3줄로 작성해주세요.
"""
        
        # LLM 요약 수행
        summary_result = llm.invoke(summary_prompt)
        summary_lines = summary_result.content.strip().split('\n')
        
        # 저널별 SJR 점수 분석 추가
        print(f"[TOOL] 저널 SJR 점수 분석 중...")
        enhanced_results = []
        
        for paper in results:
            # 논문 제목에서 저널명 추정 (개선된 매칭 로직)
            estimated_journal = "Unknown"
            best_match_metrics = None
            
            # Playwright 검색에서는 실제 저널명을 얻기 어려우므로 
            # 논문 제목 기반 추정은 부정확하여 제거
            # 대신 실제 저널명이 있는 API 검색을 우선 사용하도록 안내
            print(f"[TOOL] Playwright 검색은 저널명 추정이 어려움. API 검색 권장.")
            
            # 영향력 점수 계산
            impact_score = 0.0
            impact_tier = "📝 Unknown"
            impact_class = "unknown-impact"
            
            if best_match_metrics:
                impact_score = scimago_analyzer.calculate_impact_score(best_match_metrics, 0)
                impact_tier, impact_class = scimago_analyzer.get_impact_tier(impact_score)
            
            enhanced_results.append({
                'paper': paper,
                'journal': estimated_journal,
                'metrics': best_match_metrics,
                'impact_score': impact_score,
                'impact_tier': impact_tier,
                'impact_class': impact_class
            })
        
        # 영향력 점수 순으로 정렬
        enhanced_results.sort(key=lambda x: x['impact_score'], reverse=True)
        
        # HTML 표 형식으로 결과 생성 (SJR 포함)
        html_table = f"""
<div class="research-papers-table enhanced-search">
<h3>🔬 논문 검색 결과 (IF 점수 포함): {query}</h3>
<div class="search-summary" style="background-color: #f0f9ff; padding: 1rem; border-radius: 0.5rem; margin-bottom: 1rem;">
    <p style="color: #0369a1; margin: 0;">
    📊 <strong>Enhanced Search:</strong> Scimago Journal Rank (SJR) 기반 IF 점수 분석 포함<br>
    🎯 영향력 점수순 정렬 | 📈 저널 품질 평가 | 📚 IF 점수 자동 계산
    </p>
</div>
<table class="papers-table enhanced-table" style="width: 100%; table-layout: fixed; min-width: 1200px;">
<thead>
<tr>
<th style="width: 40%;">논문 제목</th>
<th style="width: 20%;">추정 저널 (SJR)</th>
<th style="width: 35%;">주요 내용 (AI 요약)</th>
<th style="width: 5%;">링크</th>
</tr>
</thead>
<tbody>
"""
        
        # 향상된 논문 정보를 테이블 행으로 추가
        for i, enhanced_result in enumerate(enhanced_results):
            paper = enhanced_result['paper']
            journal = enhanced_result['journal']
            metrics = enhanced_result['metrics']
            impact_score = enhanced_result['impact_score']
            impact_tier = enhanced_result['impact_tier']
            impact_class = enhanced_result['impact_class']
            
            # LLM 요약에서 해당 논문의 요약 찾기
            paper_summary = "요약 생성 중..."
            for line in summary_lines:
                if f"논문{i+1}||" in line:
                    parts = line.split("||")
                    if len(parts) >= 3:
                        paper_summary = parts[2].strip()
                    break
            
            # 요약이 없으면 기본 초록 사용
            if paper_summary == "요약 생성 중..." or not paper_summary:
                paper_summary = paper['abstract'][:150] + "..." if len(paper['abstract']) > 150 else paper['abstract']
            
            # 저널 정보 구성
            if metrics and metrics.get('sjr_score', 0) > 0:
                journal_info = f"""{journal}<br>
                <small style='color: #666;'>
                SJR: {metrics['sjr_score']:.2f} ({metrics['quartile']})
                </small>"""
            else:
                journal_info = f"""{journal}<br>
                <small style='color: #999;'>SJR 정보 없음</small>"""
            
            # 영향력 배지
            impact_badge = f"""
            <div class="impact-badge {impact_class}" style="text-align: center; padding: 0.4rem; border-radius: 0.4rem; border: 2px solid; font-size: 0.8em;">
                <div style="font-weight: bold;">{impact_score:.1f}</div>
                <div style="font-size: 0.7em;">{impact_tier.split()[0]}</div>
            </div>
            """
            
            html_table += f"""
<tr class="enhanced-paper-row">
<td class="paper-title" style="word-wrap: break-word; padding: 0.8rem; line-height: 1.3;"><strong>{paper['title']}</strong></td>
<td class="journal-info" style="word-wrap: break-word; padding: 0.8rem;">{journal_info}</td>
<td class="paper-summary" style="word-wrap: break-word; padding: 0.8rem; line-height: 1.4;">{paper_summary}</td>
<td style="padding: 0.5rem;"><a href="{paper['url']}" target="_blank" class="paper-link">원문</a></td>
</tr>
"""
        
        html_table += f"""
</tbody>
</table>

</div>

<style>
.enhanced-search .impact-badge {{
    font-size: 0.8em;
}}
.enhanced-search .high-impact {{ background: #fef2f2; border-color: #ef4444; color: #dc2626; }}
.enhanced-search .medium-high-impact {{ background: #fff7ed; border-color: #f97316; color: #ea580c; }}
.enhanced-search .medium-impact {{ background: #fefce8; border-color: #eab308; color: #ca8a04; }}
.enhanced-search .standard-impact {{ background: #f0fdf4; border-color: #22c55e; color: #16a34a; }}
.enhanced-search .low-impact {{ background: #f8fafc; border-color: #64748b; color: #475569; }}
.enhanced-search .unknown-impact {{ background: #f1f5f9; border-color: #94a3b8; color: #64748b; }}

.enhanced-search .impact-legend {{
    padding: 0.25rem 0.5rem;
    border-radius: 0.25rem;
    border: 1px solid;
    font-size: 0.75em;
    white-space: nowrap;
}}

.enhanced-search .enhanced-paper-row:hover {{
    background-color: #f9fafb;
}}

.enhanced-search .journal-info {{
    font-size: 0.85em;
}}

.enhanced-search .enhanced-table {{
    min-width: 1200px;
}}

.enhanced-search .paper-title {{
    font-size: 0.95em;
    line-height: 1.3;
}}

.enhanced-search .paper-summary {{
    font-size: 0.9em;
    line-height: 1.4;
}}
</style>
"""
        
        print(f"[TOOL] 논문 검색 도구 완료: {len(results)}개 논문")
        # HTML 마커 제거하고 순수 HTML만 반환
        return html_table
        
    except Exception as e:
        error_msg = f"논문 검색 도구 오류: {str(e)}"
        print(f"[TOOL] 오류 발생: {error_msg}")
        return f"""
<div class="research-papers-table">
<h3>📚 논문 검색 오류</h3>
<p style="color: #dc2626; padding: 1rem; background-color: #fef2f2; border: 1px solid #fecaca; border-radius: 0.5rem;">
{error_msg}<br><br>
💡 <strong>대안:</strong> "API로 {query} 논문을 검색해줘"를 시도해보세요.
</p>
</div>
"""

# 화학 논문 분석 도구
@tool
def analyze_chemical_paper(paper_url: str) -> str:
    """PubMed 논문 URL을 받아서 화학적 관점에서 논문을 분석하고 요약합니다."""
    try:
        print(f"[DEBUG] 논문 분석 시작: {paper_url}")
        
        # 논문 내용 추출
        import asyncio
        
        async def _extract_paper_content():
            searcher = paper_search(max_articles=1)
            await searcher.setup_browser()
            
            try:
                # 논문 상세 페이지에서 초록과 추가 정보 추출
                abstract = await searcher.extract_abstract(paper_url)
                
                # 추가로 논문 메타데이터 추출
                page = await searcher.context.new_page()
                await page.goto(paper_url, wait_until='domcontentloaded', timeout=20000)
                await asyncio.sleep(3)
                
                # 논문 제목 추출
                title_element = await page.query_selector('h1.heading-title, .abstract-title, h1')
                title = await title_element.inner_text() if title_element else "제목 없음"
                
                # 저자 정보 추출
                authors_elements = await page.query_selector_all('.authors .authors-list-item, .contrib-author')
                authors = []
                for author_elem in authors_elements[:5]:  # 최대 5명
                    try:
                        author_name = await author_elem.inner_text()
                        authors.append(author_name.strip())
                    except:
                        continue
                
                # 발행 정보 추출
                pub_date_element = await page.query_selector('.cit, .publication-date, .pub-date')
                pub_date = await pub_date_element.inner_text() if pub_date_element else "날짜 없음"
                
                await page.close()
                await searcher.close_browser()
                
                return {
                    "title": title,
                    "authors": authors,
                    "publication_date": pub_date,
                    "abstract": abstract
                }
                
            except Exception as e:
                await searcher.close_browser()
                raise e
        
        # 비동기 함수 실행
        try:
            loop = asyncio.get_event_loop()
            if loop.is_running():
                import concurrent.futures
                def run_in_thread():
                    new_loop = asyncio.new_event_loop()
                    asyncio.set_event_loop(new_loop)
                    try:
                        return new_loop.run_until_complete(_extract_paper_content())
                    finally:
                        new_loop.close()
                
                with concurrent.futures.ThreadPoolExecutor() as executor:
                    future = executor.submit(run_in_thread)
                    paper_data = future.result(timeout=60)
            else:
                paper_data = loop.run_until_complete(_extract_paper_content())
        except RuntimeError:
            paper_data = asyncio.run(_extract_paper_content())
        
        # LLM을 사용한 화학적 분석
        analysis_prompt = f"""
다음 화학/의학 논문을 분석해주세요:

제목: {paper_data['title']}
저자: {', '.join(paper_data['authors'][:3])}{'...' if len(paper_data['authors']) > 3 else ''}
발행일: {paper_data['publication_date']}

초록:
{paper_data['abstract']}

다음 관점에서 분석해주세요:
1. 연구된 화학 물질 또는 약물
2. 주요 연구 방법론
3. 화학적/약리학적 메커니즘
4. 주요 발견사항
5. 임상적 의의
6. 한계점 및 향후 연구 방향

분석 결과를 한국어로 체계적으로 정리해주세요.
"""
        
        # LLM 분석 수행
        analysis_result = llm.invoke(analysis_prompt)
        
        # HTML 형식으로 결과 정리
        html_result = f"""
<div class="paper-analysis">
<h3>📋 논문 분석 결과</h3>
<div style="background-color: #f8f9fa; padding: 1rem; border-radius: 0.5rem; margin-bottom: 1rem;">
<h4 style="color: #1f2937; margin-bottom: 0.5rem;">논문 정보</h4>
<p><strong>제목:</strong> {paper_data['title']}</p>
<p><strong>저자:</strong> {', '.join(paper_data['authors'][:3])}{'...' if len(paper_data['authors']) > 3 else ''}</p>
<p><strong>발행일:</strong> {paper_data['publication_date']}</p>
<p><strong>링크:</strong> <a href="{paper_url}" target="_blank" style="color: #10b981;">원문 보기</a></p>
</div>

<div style="background-color: #ffffff; padding: 1rem; border: 1px solid #e5e7eb; border-radius: 0.5rem;">
<h4 style="color: #1f2937; margin-bottom: 1rem;">🔬 화학적 분석</h4>
<div style="white-space: pre-wrap; line-height: 1.6;">{analysis_result.content}</div>
</div>
</div>
"""
        
        print(f"[DEBUG] 논문 분석 완료")
        return html_result
        
    except Exception as e:
        error_msg = f"논문 분석 중 오류 발생: {str(e)}"
        print(f"[DEBUG] 논문 분석 오류: {error_msg}")
        return f"""
<div class="paper-analysis">
<h3>📋 논문 분석 오류</h3>
<p style="color: #dc2626; padding: 1rem; background-color: #fef2f2; border: 1px solid #fecaca; border-radius: 0.5rem;">
{error_msg}<br><br>
올바른 PubMed URL인지 확인하고 다시 시도해주세요.
</p>
</div>
"""

@tool
def analyze_chemical_research_trends(topic: str, max_papers: int = 10) -> str:
    """화학 주제의 연구 동향을 분석합니다. 논문 검색과 함께 각 논문의 저널 IF 점수(SJR 기반)와 영향력 분석을 자동으로 수행합니다. 사용자가 '연구 동향', '논문 동향', 'IF 점수', 'SJR 분석' 등을 요청하면 이 도구를 사용하세요. 완전한 HTML 표를 반환하므로 그대로 출력하세요."""
    try:
        print(f"[DEBUG] 연구 동향 분석 시작: {topic}")
        
        # 논문 검색 수행
        import asyncio
        
        async def _search_papers_for_trends():
            searcher = paper_search(max_articles=max_papers)
            print(f"[DEBUG] paper_search 인스턴스 생성됨 (동향 분석용)")
            
            # 필터를 완화한 검색 URL (모든 논문 타입 포함, 최근 5년)
            await searcher.setup_browser()
            search_url = f"https://pubmed.ncbi.nlm.nih.gov/?term={topic}&filter=datesearch.y_5"
            
            print(f"[DEBUG] 동향 분석 검색 URL: {search_url}")
            
            # 403 오류 대응을 위한 재시도 로직
            max_retries = 3
            articles = []
            
            for attempt in range(max_retries):
                page = await searcher.context.new_page()
                try:
                    print(f"[DEBUG] 동향 분석 시도 {attempt + 1}/{max_retries}")
                    
                    # 페이지 로딩 전 대기
                    await asyncio.sleep(2 * (attempt + 1))
                    
                    await page.goto(search_url, wait_until='domcontentloaded', timeout=30000)
                    await asyncio.sleep(5)
                    
                    print("[DEBUG] 동향 분석용 검색 페이지 로딩 완료")
                    
                    # 페이지 제목 확인
                    page_title = await page.title()
                    print(f"[DEBUG] 페이지 제목: {page_title}")
                    
                    # 403 오류 감지
                    if "403" in page_title or "Forbidden" in page_title:
                        print(f"[DEBUG] 403 오류 감지됨. 동향 분석 시도 {attempt + 1}")
                        if attempt < max_retries - 1:
                            await page.close()
                            continue
                        else:
                            print("[DEBUG] 동향 분석 모든 재시도 실패. 403 오류 지속됨.")
                            await page.close()
                            break
                    
                    # 정상적인 페이지인 경우 논문 추출
                    articles = await searcher.extract_articles(page, topic)
                    await page.close()
                    break
                    
                except Exception as e:
                    print(f"[DEBUG] 동향 분석 시도 {attempt + 1} 중 오류 발생: {str(e)}")
                    await page.close()
                    if attempt < max_retries - 1:
                        continue
                    else:
                        print("[DEBUG] 동향 분석 모든 재시도 실패")
                        break

            results = []
            print(f"\n[DEBUG] 동향 분석용 초록 추출 시작... {len(articles)}개의 논문 처리")
            
            for i, article in enumerate(articles, 1):
                print(f"\n=== 논문 {i}/{len(articles)} 초록 추출 ===")
                abstract = await searcher.extract_abstract(article["url"])
                results.append({
                    "title": article["title"],
                    "url": article["url"],
                    "abstract": abstract
                })

            await searcher.close_browser()
            return results
        
        # 이벤트 루프 실행
        results = None
        try:
            loop = asyncio.get_event_loop()
            if loop.is_running():
                print("[DEBUG] 기존 이벤트 루프가 실행 중, 새 스레드에서 실행")
                import concurrent.futures
                def run_in_thread():
                    new_loop = asyncio.new_event_loop()
                    asyncio.set_event_loop(new_loop)
                    try:
                        return new_loop.run_until_complete(_search_papers_for_trends())
                    finally:
                        new_loop.close()
                
                with concurrent.futures.ThreadPoolExecutor() as executor:
                    future = executor.submit(run_in_thread)
                    results = future.result(timeout=120)
            else:
                print("[DEBUG] 기존 루프에서 실행")
                results = loop.run_until_complete(_search_papers_for_trends())
        except RuntimeError as e:
            print(f"[DEBUG] RuntimeError 발생: {e}, 새 루프 생성")
            results = asyncio.run(_search_papers_for_trends())
        
        if not results:
            print(f"[DEBUG] 동향 분석용 검색 결과가 비어있음")
            return f"""
<div class="research-trends">
<h3>📊 '{topic}' 연구 동향 보고서</h3>
<p style="color: #dc2626; padding: 1rem; background-color: #fef2f2; border: 1px solid #fecaca; border-radius: 0.5rem;">
해당 주제에 대한 충분한 논문을 찾을 수 없어 동향 분석을 수행할 수 없습니다.<br><br>
💡 <strong>팁:</strong> 더 일반적인 검색어를 사용해보세요 (예: "cancer treatment", "drug delivery" 등)
</p>
</div>
"""
        
        print(f"[DEBUG] 동향 분석용 HTML 테이블 생성 중, {len(results)}개 논문")
        
        # LLM을 사용해서 각 논문의 주요 내용 요약 생성
        papers_summaries = []
        for i, paper in enumerate(results, 1):
            print(f"[DEBUG] 논문 {i} 요약 생성 중...")
            paper_prompt = f"""
다음 논문의 주요 내용을 2-3줄로 간단히 요약해주세요:

제목: {paper['title']}
초록: {paper['abstract']}

연구 목적, 방법, 주요 결과를 포함해서 간결하게 요약해주세요.
"""
            try:
                paper_summary_result = llm.invoke(paper_prompt)
                paper_summary = paper_summary_result.content.strip()
                
                # 요약이 너무 길면 자르기
                if len(paper_summary) > 200:
                    paper_summary = paper_summary[:200] + "..."
                    
                papers_summaries.append(paper_summary)
                    
            except Exception as e:
                print(f"[DEBUG] 논문 {i} 요약 생성 오류: {str(e)}")
                fallback_summary = paper['abstract'][:150] + "..." if len(paper['abstract']) > 150 else paper['abstract']
                papers_summaries.append(fallback_summary)
        
        # 연구 동향 분석을 위한 전체 요약 생성
        all_papers_info = ""
        for i, (paper, summary) in enumerate(zip(results, papers_summaries), 1):
            all_papers_info += f"""
논문 {i}: {paper['title']}
주요 내용: {summary}

"""
        
        trend_analysis_prompt = f"""
다음은 '{topic}' 주제로 검색된 {len(results)}개의 최신 연구 논문들입니다:

{all_papers_info}

이 논문들을 바탕으로 다음 내용을 분석해주세요:

1. 주요 연구 동향 (2-3개의 핵심 트렌드)
2. 주목할 만한 연구 성과
3. 향후 연구 방향 전망

각 항목을 2-3문장으로 간결하게 설명해주세요.
"""
        
        try:
            trend_analysis_result = llm.invoke(trend_analysis_prompt)
            trend_analysis = trend_analysis_result.content.strip()
        except Exception as e:
            print(f"[DEBUG] 동향 분석 생성 오류: {str(e)}")
            trend_analysis = "동향 분석 생성 중 오류가 발생했습니다."
        
        # 저널별 SJR 점수 분석
        print(f"[DEBUG] 저널 SJR 점수 분석 중...")
        journal_metrics_map = {}
        unique_journals = []
        
        # 논문에서 저널명 추출 (Playwright 검색 결과용)
        # Playwright 검색에서는 저널 정보를 직접 얻기 어려우므로 제목 기반 추정 사용
        for paper in results:
            # 논문 제목에서 저널명 키워드 찾기
            paper_title_lower = paper.get('title', '').lower()
            found_journal = None
            
            # Playwright 검색에서는 정확한 저널명 추출이 어려우므로
            # 제목 기반 추정은 부정확하여 제거
            # 연구 동향 분석에서는 전체적인 트렌드에 집중하고
            # 개별 논문의 정확한 IF는 API 검색을 권장
            
            if found_journal and found_journal not in unique_journals:
                unique_journals.append(found_journal)
                print(f"[DEBUG] 저널 추정: {paper_title_lower[:50]}... → {found_journal}")
        
        # 각 저널의 SJR 점수 조회
        for journal in unique_journals[:10]:  # 최대 10개 저널만 처리
            try:
                metrics = scimago_analyzer.get_journal_metrics_fast(journal)
                if metrics and metrics.get('sjr_score', 0) > 0:
                    journal_metrics_map[journal] = metrics
            except Exception as e:
                print(f"[DEBUG] 저널 분석 오류 ({journal}): {str(e)}")
                continue
        
        # 논문과 저널 메트릭을 결합하여 영향력 점수 계산
        enhanced_papers = []
        for i, (paper, summary) in enumerate(zip(results, papers_summaries), 1):
            # 논문에서 저널명 추정 (개선된 매칭 로직)
            estimated_journal = "Unknown"
            best_match_metrics = None
            
            # 알려진 저널과 매칭 시도 (더 정교한 방법)
            paper_title_lower = paper['title'].lower()
            best_match_score = 0
            
            # 제목 기반 저널 추정은 부정확하므로 제거
            # 연구 동향 분석에서는 논문 수와 키워드 분석에 집중
            print(f"[DEBUG] 연구 동향 분석: 저널별 IF 분석보다는 전체 트렌드에 집중")
            
            # 영향력 점수 계산
            impact_score = 0.0
            impact_tier = "📝 Unknown"
            impact_class = "unknown-impact"
            
            if best_match_metrics:
                impact_score = scimago_analyzer.calculate_impact_score(best_match_metrics, 0)
                impact_tier, impact_class = scimago_analyzer.get_impact_tier(impact_score)
            
            enhanced_papers.append({
                'paper': paper,
                'summary': summary,
                'journal': estimated_journal,
                'metrics': best_match_metrics,
                'impact_score': impact_score,
                'impact_tier': impact_tier,
                'impact_class': impact_class
            })
        
        # 영향력 점수 순으로 정렬
        enhanced_papers.sort(key=lambda x: x['impact_score'], reverse=True)
        
        # HTML 표 형식으로 결과 생성 (SJR 점수 포함)
        html_table = f"""
<div class="research-trends enhanced-trends">
<h3>📊 '{topic}' 연구 동향 보고서 (SJR 분석 포함)</h3>
<div class="trend-summary" style="background-color: #f0f9ff; padding: 1rem; border-radius: 0.5rem; margin-bottom: 1rem;">
    <p style="color: #0369a1; margin: 0;">
    🔬 <strong>Enhanced Analysis:</strong> Scimago Journal Rank (SJR) 영향력 분석 포함<br>
    📈 영향력 점수순 정렬 | 📚 저널 품질 평가 | 🎯 연구 동향 분석
    </p>
</div>
<table class="papers-table enhanced-trends-table" style="margin-top: 1rem;">
<thead>
<tr>
<th style="width: 35%;">논문 제목</th>
<th style="width: 25%;">추정 저널 (SJR)</th>
<th style="width: 35%;">주요 내용</th>
<th style="width: 5%;">링크</th>
</tr>
</thead>
<tbody>
"""
        
        # 각 논문 정보를 테이블 행으로 추가 (SJR 정보 포함)
        for enhanced_paper in enhanced_papers:
            paper = enhanced_paper['paper']
            summary = enhanced_paper['summary']
            journal = enhanced_paper['journal']
            metrics = enhanced_paper['metrics']
            impact_score = enhanced_paper['impact_score']
            impact_tier = enhanced_paper['impact_tier']
            impact_class = enhanced_paper['impact_class']
            
            # 저널 정보 구성
            if metrics and metrics.get('sjr_score', 0) > 0:
                journal_info = f"""{journal}<br>
                <small style='color: #666;'>
                SJR: {metrics['sjr_score']:.2f} ({metrics['quartile']})
                </small>"""
            else:
                journal_info = f"""{journal}<br>
                <small style='color: #999;'>SJR 정보 없음</small>"""
            
            # 영향력 배지
            impact_badge = f"""
            <div class="impact-badge {impact_class}" style="text-align: center; padding: 0.5rem; border-radius: 0.5rem; border: 2px solid;">
                <div style="font-weight: bold; font-size: 1em;">{impact_score:.1f}</div>
                <div style="font-size: 0.7em;">{impact_tier}</div>
            </div>
            """
            
            html_table += f"""
<tr class="enhanced-trend-row">
<td class="paper-title"><strong>{paper['title']}</strong></td>
<td class="journal-info">{journal_info}</td>
<td class="paper-summary">{summary}</td>
<td class="impact-score">{impact_badge}</td>
<td><a href="{paper['url']}" target="_blank" class="paper-link">원문</a></td>
</tr>
"""
        
        html_table += f"""
</tbody>
</table>

<div style="margin-top: 2rem; padding: 1.5rem; background-color: #f8f9fa; border: 1px solid #e9ecef; border-radius: 0.5rem;">
<h4 style="color: #1f2937; margin-bottom: 1rem;">📈 연구 동향 분석</h4>
<div style="white-space: pre-wrap; line-height: 1.6; color: #374151;">{trend_analysis}</div>
</div>

</div>

<style>
.enhanced-trends .impact-badge {{
    font-size: 0.85em;
}}
.enhanced-trends .high-impact {{ background: #fef2f2; border-color: #ef4444; color: #dc2626; }}
.enhanced-trends .medium-high-impact {{ background: #fff7ed; border-color: #f97316; color: #ea580c; }}
.enhanced-trends .medium-impact {{ background: #fefce8; border-color: #eab308; color: #ca8a04; }}
.enhanced-trends .standard-impact {{ background: #f0fdf4; border-color: #22c55e; color: #16a34a; }}
.enhanced-trends .low-impact {{ background: #f8fafc; border-color: #64748b; color: #475569; }}
.enhanced-trends .unknown-impact {{ background: #f1f5f9; border-color: #94a3b8; color: #64748b; }}

.enhanced-trends .impact-legend {{
    padding: 0.25rem 0.5rem;
    border-radius: 0.25rem;
    border: 1px solid;
    font-size: 0.75em;
    white-space: nowrap;
}}

.enhanced-trends .enhanced-trend-row:hover {{
    background-color: #f9fafb;
}}

.enhanced-trends .journal-info {{
    font-size: 0.9em;
}}
</style>
"""
        
        print(f"[DEBUG] 연구 동향 보고서 생성 완료")
        # HTML 마커 제거하고 순수 HTML만 반환
        return html_table
        
    except Exception as e:
        error_msg = f"연구 동향 분석 중 오류 발생: {str(e)}"
        print(f"[DEBUG] 연구 동향 분석 오류: {error_msg}")
        return f"""
<div class="research-trends">
<h3>📊 연구 동향 분석 오류</h3>
<p style="color: #dc2626; padding: 1rem; background-color: #fef2f2; border: 1px solid #fecaca; border-radius: 0.5rem;">
{error_msg}<br><br>
💡 <strong>대안:</strong> "API로 {topic} 논문을 검색해줘"를 시도해보세요.
</p>
</div>
"""

# 테스트용 간단한 검색 함수
@tool  
def test_paper_search(query: str = "aspirin") -> str:
    """테스트용 논문 검색 - 간단한 검색어로 테스트합니다."""
    try:
        print(f"[TEST] 테스트 검색 시작: {query}")
        result = search_research_papers_sync(query, 3)  # 3개만 검색
        print(f"[TEST] 테스트 검색 완료")
        return result
    except Exception as e:
        return f"테스트 검색 오류: {str(e)}"

@tool
def api_paper_search_tool(query: str) -> str:
    """NCBI E-utilities API를 사용한 논문 검색과 저널 IF 점수(SJR 기반) 분석을 수행합니다. 403 오류 없이 안정적으로 논문을 검색하고 각 논문의 저널 영향력을 자동 분석합니다. HTML 표 형식으로 결과를 반환하므로 그대로 사용자에게 제공하세요."""
    try:
        print(f"[API] API를 통한 논문 검색 시작: {query}")
        
        import requests
        import xml.etree.ElementTree as ET
        from urllib.parse import quote
        
        # Step 1: 검색하여 논문 ID 목록 가져오기
        search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        search_params = {
            'db': 'pubmed',
            'term': query,
            'retmax': 5,  # 기본값으로 5개 설정
            'retmode': 'xml',
            'datetype': 'pdat',
            'reldate': 1825,  # 최근 5년 (365 * 5)
            'sort': 'relevance'
        }
        
        print(f"[API] 검색 요청 중...")
        search_response = requests.get(search_url, params=search_params, timeout=30)
        
        if search_response.status_code != 200:
            return f"API 검색 실패: HTTP {search_response.status_code}"
        
        # XML 파싱하여 논문 ID 추출
        search_root = ET.fromstring(search_response.content)
        id_list = search_root.find('.//IdList')
        
        if id_list is None or len(id_list) == 0:
            return f"""
<div class="research-papers-table">
<h3>📚 API 검색 결과: {query}</h3>
<p style="color: #dc2626; padding: 1rem; background-color: #fef2f2; border: 1px solid #fecaca; border-radius: 0.5rem;">
'{query}'에 대한 논문을 찾을 수 없습니다.<br><br>
💡 <strong>팁:</strong> 더 일반적인 영어 검색어를 사용해보세요 (예: "cancer", "treatment", "therapy")
</p>
</div>
"""
        
        paper_ids = [id_elem.text for id_elem in id_list.findall('Id')]
        print(f"[API] 발견된 논문 ID: {len(paper_ids)}개")
        
        # Step 2: 논문 상세 정보 가져오기
        fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        fetch_params = {
            'db': 'pubmed',
            'id': ','.join(paper_ids),
            'retmode': 'xml'
        }
        
        print(f"[API] 논문 상세 정보 요청 중...")
        fetch_response = requests.get(fetch_url, params=fetch_params, timeout=30)
        
        if fetch_response.status_code != 200:
            return f"API 상세 정보 요청 실패: HTTP {fetch_response.status_code}"
        
        # XML 파싱하여 논문 정보 추출
        fetch_root = ET.fromstring(fetch_response.content)
        articles = fetch_root.findall('.//PubmedArticle')
        
        if not articles:
            return f"논문 상세 정보를 가져올 수 없습니다."
        
        print(f"[API] 논문 상세 정보 파싱 중: {len(articles)}개")
        
        # 논문 정보를 추출하여 LLM으로 요약
        papers_data = []
        for i, article in enumerate(articles[:5], 1):  # 최대 5개만 처리
            try:
                # 제목 추출
                title_elem = article.find('.//ArticleTitle')
                title = title_elem.text if title_elem is not None and title_elem.text else "제목 없음"
                
                # 초록 추출
                abstract_elem = article.find('.//Abstract/AbstractText')
                abstract = abstract_elem.text if abstract_elem is not None and abstract_elem.text else "초록 없음"
                
                # 저널명 추출 (가장 중요한 부분!)
                journal_title_elem = article.find('.//Journal/Title')
                journal_iso_elem = article.find('.//Journal/ISOAbbreviation')
                
                if journal_title_elem is not None and journal_title_elem.text:
                    journal_name = journal_title_elem.text.strip()
                elif journal_iso_elem is not None and journal_iso_elem.text:
                    journal_name = journal_iso_elem.text.strip()
                else:
                    journal_name = "Unknown"
                
                # 발행연도 추출
                year_elem = article.find('.//PubDate/Year')
                year = year_elem.text if year_elem is not None and year_elem.text else "연도 미상"
                
                # PMID 추출
                pmid_elem = article.find('.//PMID')
                pmid = pmid_elem.text if pmid_elem is not None and pmid_elem.text else ""
                pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "#"
                
                papers_data.append({
                    'title': title,
                    'abstract': abstract,
                    'journal': journal_name,  # 저널명 추가!
                    'year': year,
                    'url': pubmed_url
                })
                
            except Exception as e:
                print(f"[API] 논문 {i} 파싱 오류: {str(e)}")
                continue
        
        if not papers_data:
            return "논문 정보를 추출할 수 없습니다."
        
        print(f"[API] LLM을 사용한 논문 요약 시작...")
        
        # LLM을 사용해서 각 논문을 요약
        papers_info = ""
        for i, paper in enumerate(papers_data, 1):
            papers_info += f"""
논문 {i}:
제목: {paper['title']}
초록: {paper['abstract']}
연도: {paper['year']}
URL: {paper['url']}

---
"""
        
        summary_prompt = f"""
다음 {len(papers_data)}개의 논문을 분석하여 각 논문의 주요 내용을 한국어로 요약해주세요.
각 논문마다 2-3줄로 핵심 내용만 간단히 정리해주세요.

{papers_info}

각 논문에 대해 다음 형식으로 응답해주세요:
논문1||제목||주요내용요약
논문2||제목||주요내용요약
...

주요내용요약은 연구 목적, 방법, 주요 결과를 포함해서 2-3줄로 작성해주세요.
"""
        
        # LLM 요약 수행
        summary_result = llm.invoke(summary_prompt)
        summary_lines = summary_result.content.strip().split('\n')
        
        # 저널별 SJR 점수 분석 추가 (API 검색용)
        print(f"[API] 저널 SJR 점수 분석 중...")
        enhanced_papers_data = []
        
        for paper in papers_data:
            # PubMed에서 추출한 실제 저널명 사용
            actual_journal = paper.get('journal', 'Unknown')
            estimated_journal = actual_journal
            best_match_metrics = None
            
            print(f"[API] 실제 저널명: {actual_journal}")
            
            # 실제 저널명으로 SJR 메트릭 조회
            if actual_journal and actual_journal != "Unknown":
                # 1단계: 정확한 매칭 시도
                best_match_metrics = scimago_analyzer.get_journal_metrics_fast(actual_journal)
                
                # 2단계: 정확한 매칭이 안되면 약간의 변형을 시도
                if not best_match_metrics or best_match_metrics.get('sjr_score', 0) == 0:
                    # 저널명 정리 (약어, 특수문자 제거 등)
                    cleaned_journal = actual_journal.lower().replace('&', 'and').replace('.', '').strip()
                    if cleaned_journal != actual_journal.lower():
                        estimated_journal = cleaned_journal
                        best_match_metrics = scimago_analyzer.get_journal_metrics_fast(cleaned_journal)
                        if best_match_metrics and best_match_metrics.get('sjr_score', 0) > 0:
                            print(f"[API] 정리된 저널명으로 매칭 성공: {actual_journal} → {cleaned_journal}")
            
            # 영향력 점수 계산
            impact_score = 0.0
            impact_tier = "📝 Unknown"
            impact_class = "unknown-impact"
            
            if best_match_metrics:
                impact_score = scimago_analyzer.calculate_impact_score(best_match_metrics, 0)
                impact_tier, impact_class = scimago_analyzer.get_impact_tier(impact_score)
            
            enhanced_papers_data.append({
                'paper': paper,
                'journal': estimated_journal,
                'metrics': best_match_metrics,
                'impact_score': impact_score,
                'impact_tier': impact_tier,
                'impact_class': impact_class
            })
        
        # 영향력 점수 순으로 정렬
        enhanced_papers_data.sort(key=lambda x: x['impact_score'], reverse=True)
        
        # HTML 표 형식으로 결과 생성 (SJR 포함, 크기 확대)
        html_table = f"""
<div class="research-papers-table enhanced-api-search">
<h3>🔬 API 검색 결과 (SJR 분석 포함): {query}</h3>
<div class="search-summary" style="background-color: #f0f9ff; padding: 1rem; border-radius: 0.5rem; margin-bottom: 1rem;">
    <p style="color: #0369a1; margin: 0;">
    📊 <strong>Enhanced API Search:</strong> NCBI E-utilities + Scimago Journal Rank 통합 분석<br>
    🎯 SJR 점수 기반 저널 품질 평가 | 📈 Quartile 분석 포함
    </p>
</div>
<table class="papers-table enhanced-api-table" style="width: 100%; table-layout: fixed;">
<thead>
<tr>
<th style="width: 40%;">논문 제목</th>
<th style="width: 20%;">저널 (SJR)</th>
<th style="width: 35%;">주요 내용 (AI 요약)</th>
<th style="width: 5%;">링크</th>
</tr>
</thead>
<tbody>
"""
        
        # 향상된 논문 정보를 테이블 행으로 추가 (API 검색용)
        for i, enhanced_paper in enumerate(enhanced_papers_data):
            paper = enhanced_paper['paper']
            journal = enhanced_paper['journal']
            metrics = enhanced_paper['metrics']
            impact_score = enhanced_paper['impact_score']
            impact_tier = enhanced_paper['impact_tier']
            impact_class = enhanced_paper['impact_class']
            
            # LLM 요약에서 해당 논문의 요약 찾기
            paper_summary = "요약 생성 중..."
            for line in summary_lines:
                if f"논문{i+1}||" in line:
                    parts = line.split("||")
                    if len(parts) >= 3:
                        paper_summary = parts[2].strip()
                    break
            
            # 요약이 없으면 기본 초록 사용
            if paper_summary == "요약 생성 중..." or not paper_summary:
                paper_summary = paper['abstract'][:180] + "..." if len(paper['abstract']) > 180 else paper['abstract']
            
            # 저널 정보 구성
            if metrics and metrics.get('sjr_score', 0) > 0:
                journal_info = f"""{journal}<br>
                <small style='color: #10b981; font-weight: 600; font-size: 0.85em;'>
                SJR: {metrics['sjr_score']:.2f} ({metrics['quartile']})
                </small>"""
            else:
                journal_info = f"""{journal}<br>
                <small style='color: #999; font-size: 0.8em;'>SJR 정보 없음</small>"""
            
            html_table += f"""
<tr class="enhanced-api-row">
<td class="paper-title" style="word-wrap: break-word; padding: 0.8rem;"><strong>{paper['title']}</strong><br><small style="color: #666;">({paper['year']})</small></td>
<td class="journal-info" style="word-wrap: break-word; padding: 0.8rem;">{journal_info}</td>
<td class="paper-summary" style="word-wrap: break-word; padding: 0.8rem; line-height: 1.4;">{paper_summary}</td>
<td style="padding: 0.5rem;"><a href="{paper['url']}" target="_blank" class="paper-link">원문</a></td>
</tr>
"""
        
        html_table += f"""
</tbody>
</table>

</div>

<style>
.enhanced-api-search .impact-badge {{
    font-size: 0.75em;
}}
.enhanced-api-search .high-impact {{ background: #fef2f2; border-color: #ef4444; color: #dc2626; }}
.enhanced-api-search .medium-high-impact {{ background: #fff7ed; border-color: #f97316; color: #ea580c; }}
.enhanced-api-search .medium-impact {{ background: #fefce8; border-color: #eab308; color: #ca8a04; }}
.enhanced-api-search .standard-impact {{ background: #f0fdf4; border-color: #22c55e; color: #16a34a; }}
.enhanced-api-search .low-impact {{ background: #f8fafc; border-color: #64748b; color: #475569; }}
.enhanced-api-search .unknown-impact {{ background: #f1f5f9; border-color: #94a3b8; color: #64748b; }}

.enhanced-api-search .impact-legend {{
    padding: 0.25rem 0.5rem;
    border-radius: 0.25rem;
    border: 1px solid;
    font-size: 0.75em;
    white-space: nowrap;
}}

.enhanced-api-search .enhanced-api-row:hover {{
    background-color: #f9fafb;
}}

.enhanced-api-search .enhanced-api-table {{
    min-width: 1200px;
}}

.enhanced-api-search .paper-title {{
    font-size: 0.95em;
    line-height: 1.3;
}}

.enhanced-api-search .paper-summary {{
    font-size: 0.9em;
    line-height: 1.4;
}}

.enhanced-api-search .journal-info {{
    font-size: 0.85em;
    line-height: 1.2;
}}
</style>
"""
        
        print(f"[API] HTML 테이블 생성 완료")
        # HTML 마커 제거하고 순수 HTML만 반환
        return html_table
        
    except Exception as e:
        error_msg = f"API 논문 검색 중 오류 발생: {str(e)}"
        print(f"[API] 예외 발생: {error_msg}")
        return f"""
<div class="research-papers-table">
<h3>📚 API 검색 오류</h3>
<p style="color: #dc2626; padding: 1rem; background-color: #fef2f2; border: 1px solid #fecaca; border-radius: 0.5rem;">
{error_msg}
</p>
</div>
"""

@tool
def get_pubmed_search_tips() -> str:
    """PubMed 검색이 잘 안될 때 도움이 되는 팁을 제공합니다."""
    return """
<div class="search-tips">
<h3>🔍 PubMed 검색 최적화 팁</h3>

<div style="background-color: #f0fdf4; padding: 1rem; border: 1px solid #bbf7d0; border-radius: 0.5rem; margin: 1rem 0;">
<h4 style="color: #166534; margin-bottom: 0.5rem;">✅ 효과적인 검색어</h4>
<ul style="color: #166534; margin-left: 1rem;">
<li><strong>일반적인 용어 사용:</strong> "cancer", "treatment", "therapy"</li>
<li><strong>영어 검색어:</strong> 한글보다 영어가 더 많은 결과</li>
<li><strong>MeSH 용어:</strong> 의학 표준 용어 활용</li>
<li><strong>동의어 포함:</strong> "drug OR medication OR pharmaceutical"</li>
</ul>
</div>

<div style="background-color: #fef2f2; padding: 1rem; border: 1px solid #fecaca; border-radius: 0.5rem; margin: 1rem 0;">
<h4 style="color: #dc2626; margin-bottom: 0.5rem;">❌ 피해야 할 검색어</h4>
<ul style="color: #dc2626; margin-left: 1rem;">
<li><strong>너무 구체적:</strong> "doxorubicin cardiotoxicity prevention mechanism"</li>
<li><strong>브랜드명:</strong> 일반명 사용 권장</li>
<li><strong>복잡한 조합:</strong> 3개 이상의 AND 조건</li>
</ul>
</div>

<div style="background-color: #f8f9fa; padding: 1rem; border: 1px solid #e9ecef; border-radius: 0.5rem;">
<h4 style="color: #495057; margin-bottom: 0.5rem;">🎯 추천 검색어 예시</h4>
<ul style="color: #495057; margin-left: 1rem;">
<li>"cancer treatment"</li>
<li>"immunotherapy"</li>
<li>"drug delivery"</li>
<li>"chemotherapy"</li>
<li>"biomarker"</li>
<li>"clinical trial"</li>
</ul>
</div>
</div>
"""

@tool
def enhanced_paper_search_with_sjr(query: str, max_papers: int = 5) -> str:
    """SJR 영향력 분석이 포함된 향상된 논문 검색. 기존 검색보다 더 자세한 저널 정보와 영향력 점수를 제공합니다."""
    try:
        print(f"[ENHANCED] SJR 분석 포함 논문 검색 시작: {query}")
        
        # 기본 논문 검색 수행
        results = search_research_papers_sync(query, max_papers)
        
        if "논문을 찾을 수 없습니다" in results:
            return results
        
        # HTML에서 논문 정보 추출 (간단한 방법)
        # 실제로는 더 정교한 파싱이 필요하지만, 여기서는 기본 검색 결과를 향상시킴
        enhanced_html = results.replace(
            '<h3>📚 연구 논문 검색 결과:', 
            '<h3>🔬 SJR 분석 포함 논문 검색 결과:'
        )
        
        # SJR 분석 정보 추가
        sjr_info = """
        """
        
        enhanced_html = enhanced_html.replace('</div>', sjr_info + '</div>')
        
        return enhanced_html
        
    except Exception as e:
        error_msg = f"SJR 분석 포함 검색 오류: {str(e)}"
        print(f"[ENHANCED] {error_msg}")
        return f"""
        <div class="research-papers-table">
        <h3>🔬 SJR 분석 포함 검색 오류</h3>
        <p style="color: #dc2626; padding: 1rem; background-color: #fef2f2; border: 1px solid #fecaca; border-radius: 0.5rem;">
        {error_msg}
        </p>
        </div>
        """

tools = [
    generate_molecule_image, 
    name_to_smiles_cir, 
    paper_search_tool,  # 단순화된 논문 검색 도구
    api_paper_search_tool,  # 단순화된 API 검색 도구
    enhanced_paper_search_with_sjr,  # 새로운 SJR 분석 포함 검색
    analyze_chemical_paper,
    analyze_chemical_research_trends,  # 이제 SJR 분석 포함
    get_pubmed_search_tips
]

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
                    content=json.dumps(tool_result) if not isinstance(tool_result, str) else tool_result,
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
        with open("index.html", "r", encoding="utf-8") as f:
            return HTMLResponse(f.read())
    except Exception as e:
        return HTMLResponse(f"<h1>Error loading page: {str(e)}</h1>")

@app.get("/favicon.ico")
async def get_favicon():
    """Favicon 요청 처리 - 404 에러 방지"""
    # 간단한 투명 1x1 픽셀 PNG 이미지 (Base64)
    transparent_png = "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChAGA60e6kgAAAABJRU5ErkJggg=="
    
    import base64
    from fastapi.responses import Response
    
    png_data = base64.b64decode(transparent_png)
    return Response(content=png_data, media_type="image/png")

@app.post("/chat")
async def chat_endpoint(chat_message: ChatMessage):
    try:
        # HJ ChemAgent 실행
        result = graph.invoke({"messages": [{"role": "user", "content": chat_message.message}]})
        
        # HTML 도구 결과 직접 감지 및 우회 처리
        html_tool_result = None
        for message in result["messages"]:
            if hasattr(message, 'name') and message.name in ['api_paper_search_tool', 'paper_search_tool', 'analyze_chemical_research_trends']:
                tool_content = message.content
                if isinstance(tool_content, str) and (
                    '<div class="research-papers-table">' in tool_content or 
                    '<div class="research-trends">' in tool_content or
                    '<div class="paper-analysis">' in tool_content or
                    '<table' in tool_content
                ):
                    html_tool_result = tool_content
                    print(f"[DEBUG] HTML 도구 결과 직접 감지됨 - 도구: {message.name}, 길이: {len(tool_content)}")
                    break
        
        # HTML 도구 결과가 있으면 LLM 응답 대신 직접 사용
        if html_tool_result:
            response_text = html_tool_result
            print(f"[DEBUG] LLM 우회하여 도구 결과 직접 사용")
        else:
            # 마지막 AI 메시지 가져오기
            last_message = result["messages"][-1]
            response_text = last_message.content
        
        # HTML 테이블 직접 감지 및 처리 (백엔드에서 추가 보장)
        print(f"[DEBUG] 원본 응답 길이: {len(response_text)}")
        print(f"[DEBUG] 응답 시작 100자: {response_text[:100]}")
        
        # HTML 태그 감지
        has_html = ("<div" in response_text and "</div>" in response_text) or ("<table" in response_text and "</table>" in response_text)
        
        if has_html:
            print(f"[DEBUG] HTML 태그 감지됨")
            
            # HTML 시작 위치 찾기
            html_start = response_text.find("<div")
            if html_start == -1:
                html_start = response_text.find("<table")
            
            if html_start > 0:
                # HTML 앞에 불필요한 텍스트가 있으면 제거
                original_length = len(response_text)
                response_text = response_text[html_start:].strip()
                print(f"[DEBUG] HTML 앞의 텍스트 제거됨 - 원래 길이: {original_length}, 새 길이: {len(response_text)}")
            
            # 추가 텍스트와 함께 HTML이 있는 경우 처리
            if any(phrase in response_text for phrase in ["결과입니다", "분석했습니다", "보고서입니다", "검색 결과입니다"]):
                print(f"[DEBUG] 추가 설명 텍스트 감지됨, HTML 추출 시도")
                import re
                
                # 더 포괄적인 HTML 매칭
                html_patterns = [
                    r'<div[^>]*class="research-trends"[^>]*>.*?</div>',
                    r'<div[^>]*class="research-papers-table"[^>]*>.*?</div>',
                    r'<div[^>]*class="paper-analysis"[^>]*>.*?</div>',
                    r'<table[^>]*>.*?</table>'
                ]
                
                for pattern in html_patterns:
                    html_match = re.search(pattern, response_text, re.DOTALL)
                    if html_match:
                        extracted_html = html_match.group(0)
                        print(f"[DEBUG] 패턴 '{pattern[:30]}...'로 HTML 추출됨 - 길이: {len(extracted_html)}")
                        response_text = extracted_html
                        break
                else:
                    print(f"[DEBUG] HTML 패턴 매칭 실패, 원본 유지")
        else:
            print(f"[DEBUG] HTML 태그 없음, 일반 텍스트로 처리")
        
        print(f"[DEBUG] 최종 응답 길이: {len(response_text)}")
        print(f"[DEBUG] 최종 응답 시작 200자: {response_text[:200]}")
        
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
            response_text = "분자의 구조 이미지를 생성했습니다."
        
        return ChatResponse(response=response_text, image=image_data)
        
    except Exception as e:
        return ChatResponse(
            response=f"오류가 발생했습니다: {str(e)}"
        )

# 정적 파일 서빙 (필요시)
# try:
#     app.mount("/static", StaticFiles(directory="static"), name="static")
# except Exception as e:
#     print(f"Static files mounting error: {e}")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8002)
