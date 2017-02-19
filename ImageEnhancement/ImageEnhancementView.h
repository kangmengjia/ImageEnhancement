
// ImageEnhancementView.h : CImageEnhancementView ��Ľӿ�
//

#pragma once


class CImageEnhancementView : public CView
{
protected: // �������л�����
	CImageEnhancementView();
	DECLARE_DYNCREATE(CImageEnhancementView)

// ����
public:
	CImageEnhancementDoc* GetDocument() const;

	CString m_img_path;

private:
	cv::Mat m_srcImg;
	cv::Mat m_dstImg;

// ����
public:

// ��д
public:
	virtual void OnDraw(CDC* pDC);  // ��д�Ի��Ƹ���ͼ
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
protected:
	virtual BOOL OnPreparePrinting(CPrintInfo* pInfo);
	virtual void OnBeginPrinting(CDC* pDC, CPrintInfo* pInfo);
	virtual void OnEndPrinting(CDC* pDC, CPrintInfo* pInfo);

// ʵ��
public:
	virtual ~CImageEnhancementView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// ���ɵ���Ϣӳ�亯��
protected:
	afx_msg void OnFilePrintPreview();
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnContextMenu(CWnd* pWnd, CPoint point);
	DECLARE_MESSAGE_MAP()
public:
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
	afx_msg void OnFileOpen();
	afx_msg void OnRetinex();
};

#ifndef _DEBUG  // ImageEnhancementView.cpp �еĵ��԰汾
inline CImageEnhancementDoc* CImageEnhancementView::GetDocument() const
   { return reinterpret_cast<CImageEnhancementDoc*>(m_pDocument); }
#endif

