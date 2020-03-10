
// OpenCVMFCDlg.h: 头文件
//

#pragma once

#include "CEditEx.h"

#include "ImageProcess.h"
#include "Resource.h"

// COpenCVMFCDlg 对话框
class COpenCVMFCDlg : public CDialogEx {
private:
	const int MAX_ROTATE_ANGLE = 360;
	const double MAX_SCLAE_RATE = 5;
	const double MIN_SCLAE_RATE = 0.01;

	const int MAX_THREAD_COUNT = 32;

	static const int UPDATE_INTERVAL = 100;

	// 构造
public:
	COpenCVMFCDlg(CWnd* pParent = NULL);	// 标准构造函数
	~COpenCVMFCDlg();

	// 对话框数据
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_OPENCVMFC_DIALOG };
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX) {
		CDialogEx::DoDataExchange(pDX);
		DDX_Control(pDX, PATH_INPUT, m_pathEdit);
		DDX_Control(pDX, SRC_IMAGE, m_srcImage);
		DDX_Control(pDX, OUT_IMAGE, m_outImage);
		DDX_Control(pDX, ROTATE_SLIDER, m_rotateSlider);
		DDX_Control(pDX, XSCALE_SILDER, m_xSlider);
		DDX_Control(pDX, YSCALE_SILDER, m_ySlider);
		DDX_Control(pDX, ROTATE_INPUT, m_rotateInput);
		DDX_Control(pDX, XSCALE_INPUT, m_xInput);
		DDX_Control(pDX, YSCALE_INPUT, m_yInput);
		DDX_Control(pDX, SCALE_ALGO, m_scaleAlgo);
		DDX_Control(pDX, THREAD_COUNT, m_threadInput);
		DDX_Control(pDX, THREAD_MODE, m_threadMode);
		DDX_Control(pDX, DFT_ENABLE, m_dftEnable);
		DDX_Control(pDX, GAUS_AVG, m_gausAvg);
		DDX_Control(pDX, GAUS_STD, m_gausStd);
		DDX_Control(pDX, NO_NOISE, m_noNoise);
		DDX_Control(pDX, GAUS_NOISE, m_gausNoise);
		DDX_Control(pDX, SALT_NOISE, m_saltNoise);
		DDX_Control(pDX, SALT_FACTOR, m_saltFactor);
		DDX_Control(pDX, AUTO_FIT, m_autoFit);
		DDX_Control(pDX, FILTER_TYPE, m_filterType);
		DDX_Control(pDX, AUTO_REFRESH, m_autoRefresh);
		DDX_Control(pDX, GAUS_FILTER_STD, m_gausFilterStd);
	}	// DDX/DDV 支持

private:
	ProcessParam param;

	bool initializing = true;
	bool processing = false;
	bool dirty = true;

	bool autoFit = false;
	bool autoRefresh = false;

	static UINT updateDialog(void* p);

	void initialize();

	void initializeScaleAlgo();
	void initializeThreadInputs();
	void initializeFilterInputs();

	void updateDialog();

	void openConsole();
	void closeConsole();

	void refreshTransformInputs();
	void refreshNoiseInputs();

	void refreshSrcPicture();
	void refreshOutPicture();
	void refreshTransfrom();

	void refreshTransfromSingleThread();
	void refreshTransfromMultWINThread(int count);
	void refreshTransfromMultOpenMPThread(int count);

	LRESULT onProcessTranslateCompleted(WPARAM wParam, LPARAM lParam);

	// 实现
protected:

	HICON m_hIcon;
	CImage* m_pImgSrc;
	CImage* m_pImgDist;

	// 生成的消息映射函数
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()

public:
	// 图片路径
	CEdit m_pathEdit;
	CStatic m_srcImage;
	CStatic m_outImage;
	CSliderCtrl m_rotateSlider;
	CSliderCtrl m_xSlider;
	CSliderCtrl m_ySlider;
	CEditEx m_rotateInput;
	CEditEx m_xInput;
	CEditEx m_yInput;
	CComboBox m_scaleAlgo;
	CComboBox m_threadInput;
	CComboBox m_threadMode;

	afx_msg void OnBnClickedFile();

	afx_msg void OnRotateSliderChanged(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnXScaleSliderChanged(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnYScaleSliderChanged(NMHDR *pNMHDR, LRESULT *pResult);

	afx_msg void OnRotateInputBlur();
	afx_msg void OnXScaleInputBlur();
	afx_msg void OnYScaleInputBlur();

	afx_msg void OnCbnSelchangeAlgo();
	CButton m_dftEnable;
	afx_msg void OnBnClickedEnable();
	CEditEx m_gausAvg;
	CEditEx m_gausStd;
	CEditEx m_saltFactor;

	CButton m_noNoise;
	CButton m_gausNoise;
	CButton m_saltNoise;
	afx_msg void OnNoNoise();
	afx_msg void OnGausNoise();
	afx_msg void OnSaltNoise();
	afx_msg void OnEnKillfocusAvg();
	afx_msg void OnEnKillfocusStd();
	afx_msg void OnEnKillfocusFactor();
	CButton m_autoFit;
	afx_msg void OnBnClickedFit();
	CComboBox m_filterType;
	afx_msg void OnCbnSelchangeType();
	afx_msg void OnRefreshButtonClicked();
	CButton m_autoRefresh;
	afx_msg void OnBnClickedRefresh();
	CEdit m_gausFilterStd;
	afx_msg void OnEnKillfocusFilterStd();
	afx_msg void OnCbnSelchangeMode();
};
