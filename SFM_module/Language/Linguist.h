#pragma once

#include <QObject>
#include <qtranslator.h>

enum Language
{
    zh_cn,
    en_us,
    zhc_cn
};

class CLinguist :public QObject
{
    Q_OBJECT
private:
    CLinguist();

public:
    // 当前语言 简体中文、英文、未定义
    Language m_CurrentLanguage = Language::en_us;
    // 切换语言
    void ChangeLanguage(Language lan);
    static  CLinguist* GetLinguistInstance();
private:
    static CLinguist* linguist;
    QTranslator* m_trans;
signals:
    // 语言切换信号，通知当前系统语言已经被切换
    void LanguageChaned();
};
