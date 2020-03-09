#include "stdafx.h"
#include "CEditEx.h"

// CEditEx

IMPLEMENT_DYNAMIC(CEditEx, CEdit)

CEditEx::CEditEx() {}

CEditEx::~CEditEx() {}

BEGIN_MESSAGE_MAP(CEditEx, CEdit)
	ON_WM_CHAR()
	ON_WM_KEYUP()
END_MESSAGE_MAP()

// CEditEx ��Ϣ�������
void CEditEx::OnChar(UINT nChar, UINT nRepCnt, UINT nFlags) {
	// TODO:  �ڴ������Ϣ�����������/�����Ĭ��ֵ

	// ����С����
	if (nChar == '.') {
		CString str;
		GetWindowText(str);

		// ���Ƶ�һλΪС��
		if (str.GetLength() == 0) {
			// ��һλ����С����
			MessageBox(_T("��һλ��������С����"));
			return;
		}
		// ����ֻ������һ��С����
		if (str.Find('.') == -1) {
			CEdit::OnChar(nChar, nRepCnt, nFlags);
		} else {
			CString str;
			GetWindowText(str);
			if (str[0] == '.') {
				SetWindowText(str.Mid(1, str.GetLength()));
				MessageBox(_T("��һλ��������С����"));
			}
			// С������ֵڶ���
			MessageBox(_T("С����ֻ������һ��"));
		}
	}
	// �������ֺ��˸��
	else if ((nChar >= '0' && nChar <= '9') || nChar == 0x08) {
		CEdit::OnChar(nChar, nRepCnt, nFlags);
	} else {
		// ���ַ����ּ����˸��
		MessageBox(_T("ֻ���������֣��˸��"));
	}
}

// �޸�����������֮�󣬿����ڵ�һλ����С����
void CEditEx::OnKeyUp(UINT nChar, UINT nRepCnt, UINT nFlags) {
	// TODO: �ڴ������Ϣ�����������/�����Ĭ��ֵ
	if (nChar == VK_DECIMAL || nChar == VK_OEM_PERIOD) {
		CString str;
		GetWindowText(str);
		if (str[0] == '.') {
			SetWindowText(str.Mid(1, str.GetLength()));
			MessageBox(_T("��һλ��������С����"));
		}
	}
	CEdit::OnKeyUp(nChar, nRepCnt, nFlags);

}
