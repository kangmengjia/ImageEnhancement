
// ImageEnhancementDoc.cpp : CImageEnhancementDoc ���ʵ��
//

#include "stdafx.h"
// SHARED_HANDLERS ������ʵ��Ԥ��������ͼ������ɸѡ�������
// ATL ��Ŀ�н��ж��壬�����������Ŀ�����ĵ����롣
#ifndef SHARED_HANDLERS
#include "ImageEnhancement.h"
#endif

#include "ImageEnhancementDoc.h"

#include <propkey.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CImageEnhancementDoc

IMPLEMENT_DYNCREATE(CImageEnhancementDoc, CDocument)

BEGIN_MESSAGE_MAP(CImageEnhancementDoc, CDocument)
END_MESSAGE_MAP()


// CImageEnhancementDoc ����/����

CImageEnhancementDoc::CImageEnhancementDoc()
{
	// TODO:  �ڴ����һ���Թ������

}

CImageEnhancementDoc::~CImageEnhancementDoc()
{
}

BOOL CImageEnhancementDoc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;

	// TODO:  �ڴ�������³�ʼ������
	// (SDI �ĵ������ø��ĵ�)

	return TRUE;
}




// CImageEnhancementDoc ���л�

void CImageEnhancementDoc::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO:  �ڴ���Ӵ洢����
	}
	else
	{
		// TODO:  �ڴ���Ӽ��ش���
	}
}

#ifdef SHARED_HANDLERS

// ����ͼ��֧��
void CImageEnhancementDoc::OnDrawThumbnail(CDC& dc, LPRECT lprcBounds)
{
	// �޸Ĵ˴����Ի����ĵ�����
	dc.FillSolidRect(lprcBounds, RGB(255, 255, 255));

	CString strText = _T("TODO: implement thumbnail drawing here");
	LOGFONT lf;

	CFont* pDefaultGUIFont = CFont::FromHandle((HFONT) GetStockObject(DEFAULT_GUI_FONT));
	pDefaultGUIFont->GetLogFont(&lf);
	lf.lfHeight = 36;

	CFont fontDraw;
	fontDraw.CreateFontIndirect(&lf);

	CFont* pOldFont = dc.SelectObject(&fontDraw);
	dc.DrawText(strText, lprcBounds, DT_CENTER | DT_WORDBREAK);
	dc.SelectObject(pOldFont);
}

// ������������֧��
void CImageEnhancementDoc::InitializeSearchContent()
{
	CString strSearchContent;
	// ���ĵ����������������ݡ�
	// ���ݲ���Ӧ�ɡ�;���ָ�

	// ����:     strSearchContent = _T("point;rectangle;circle;ole object;")��
	SetSearchContent(strSearchContent);
}

void CImageEnhancementDoc::SetSearchContent(const CString& value)
{
	if (value.IsEmpty())
	{
		RemoveChunk(PKEY_Search_Contents.fmtid, PKEY_Search_Contents.pid);
	}
	else
	{
		CMFCFilterChunkValueImpl *pChunk = NULL;
		ATLTRY(pChunk = new CMFCFilterChunkValueImpl);
		if (pChunk != NULL)
		{
			pChunk->SetTextValue(PKEY_Search_Contents, value, CHUNK_TEXT);
			SetChunkValue(pChunk);
		}
	}
}

#endif // SHARED_HANDLERS

// CImageEnhancementDoc ���

#ifdef _DEBUG
void CImageEnhancementDoc::AssertValid() const
{
	CDocument::AssertValid();
}

void CImageEnhancementDoc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG


// CImageEnhancementDoc ����


//BOOL CImageEnhancementDoc::OnOpenDocument(LPCTSTR lpszPathName)
//{
//	if (!CDocument::OnOpenDocument(lpszPathName))
//		return FALSE;
//
//	m_Img = cv::imread((char*)lpszPathName);
//
//	DeleteContents();//ɾ���ĵ������ݣ�ȷ������ͼ������֮ǰ�ĵ�Ϊ��
//
//	SetPathName(lpszPathName);//�����ļ�����
//	SetModifiedFlag(FALSE);
//
//	return TRUE;
//}


cv::Mat CImageEnhancementDoc::GetImage()
{
	return m_Img;
}