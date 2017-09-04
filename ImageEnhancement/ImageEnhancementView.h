
// ImageEnhancementView.h : CImageEnhancementView ��Ľӿ�
//
#ifndef __ImageEnhancementView_h__
#define __ImageEnhancementView_h__


#include "ImageEnhancementDoc.h"

extern void ImageProcessByDarkWithBilateralFilter(cv::Mat src, cv::Mat dst, GuideParames para, DWORD &s, DWORD &e);
extern void ImageProcessByDark(cv::Mat src, cv::Mat dst, GuideParames para);
extern void retinex(cv::Mat src, cv::Mat im_dst, RetinexParams retinex_param);
extern void retinex_GPU(cv::Mat src, cv::Mat im_dst, RetinexParams retinex_param, DWORD &s, DWORD &e);
extern void initCUDAMem(int width, int height, int channels);
extern void initTexMem(int width, int height, uint *hImage);
extern void freeTexMem();
extern void freeCUDAMem();
extern void initCUDAMem_bil(int width, int height);
extern void freeCUDAMem_bil();
extern void updateGaussian_bil(float sigma_color, float sigma_space, int radius);
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
	int m_time_use = 0;

	RetinexParams  m_retinex_param;
	GuideParames   m_guide_param;
public:

	CMFCPropertyGridProperty* m_property_retinex;
	CMFCPropertyGridProperty* m_property_anyuanse;
	CMFCPropertyGridProperty* m_property_shuangbian;

	CString getScalesMode(int mode);
	int setScalesMode(CString mode);
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
	afx_msg LRESULT OnReloadMsg(WPARAM wparam, LPARAM lparam);
protected:
public:
	afx_msg void OnDark();
	afx_msg void OnDarkWithBilateralFilter();
};

#ifndef _DEBUG  // ImageEnhancementView.cpp �еĵ��԰汾
inline CImageEnhancementDoc* CImageEnhancementView::GetDocument() const
   { return reinterpret_cast<CImageEnhancementDoc*>(m_pDocument); }
#endif

#endif // __ImageEnhancementView_h__
