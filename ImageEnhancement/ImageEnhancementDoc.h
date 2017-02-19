
// ImageEnhancementDoc.h : CImageEnhancementDoc ��Ľӿ�
//


#pragma once


class CImageEnhancementDoc : public CDocument
{
protected: // �������л�����
	CImageEnhancementDoc();
	DECLARE_DYNCREATE(CImageEnhancementDoc)

// ����
public:
	CString m_szPathName;     // ͼƬ�ļ�·��
private:
	cv::Mat m_Img;
	


// ����
public:
	cv::Mat GetImage();


// ��д
public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);
#ifdef SHARED_HANDLERS
	virtual void InitializeSearchContent();
	virtual void OnDrawThumbnail(CDC& dc, LPRECT lprcBounds);
#endif // SHARED_HANDLERS

// ʵ��
public:
	virtual ~CImageEnhancementDoc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// ���ɵ���Ϣӳ�亯��
protected:
	DECLARE_MESSAGE_MAP()

#ifdef SHARED_HANDLERS
	// ����Ϊ�����������������������ݵ� Helper ����
	void SetSearchContent(const CString& value);
#endif // SHARED_HANDLERS
public:
//	virtual BOOL OnOpenDocument(LPCTSTR lpszPathName);
};