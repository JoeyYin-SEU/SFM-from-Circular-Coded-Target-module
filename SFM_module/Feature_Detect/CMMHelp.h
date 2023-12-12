#ifndef CMMHELP_H
#define CMMHELP_H

#include <QtGui>
#include <QApplication>

enum ScaleBarType
{
	CodeScaleBar = 0,
	UnCodeScaleBar = 1,
};
enum UnitType
{
	//国际基本单位：米（m），千克（kg) ,秒（s），卡尔文（K），安培（A），摩尔（mol） ，坎德拉（cd）

	//长度单位 nm，μm，mm，m，km，
	Nanometer = 0,
	Micrometer = 1,
	Millimeter = 2,
	Meter = 3,
	Kilometer = 4,

	//质量单位 mg , g, kg
	Milligram = 5,
	Gram = 6,
	Kilogram = 7,

	//时间单位 ms，s，min，h
	Millisecond = 8,
	Second = 9,
	Minute = 10,
	Hour = 11,

	//温度单位 ℃，℉ （F=32+1.8×C），K （K = C+273.15）
	TempertureC = 12,
	TempertureF = 13,
	TempertureK = 14,

	//电流安培 nA，μA，mA，A，kA
	Nanoampere = 15,
	Microampere = 16,
	Milliampere = 17,
	Ampere = 18,
	Kiloampere = 19,

	//物质的量 mol
	Mole = 20,

	//发光强度 cd
	Candela = 21,

	//应变 με
	Microstrain = 22,

	//无单位
	NoUnit = 23,

	//体积单位
	CUMM = 24,
	CUCM = 25,
	CUDM = 26,
	CUM = 27,
};
enum MarkPointColorType
{
	BlackDownWhiteUp = 0,
	WhiteDownBlackUp = 1,
	Uncertainty = 2,
};
enum CodePointBitesType
{
	CodeBites15 = 0,//现在是15位的
	CodeBites12 = 1,
};
enum CMMMode
{
	CalculationMode = 0,
	EvaluationMode = 1,
};
enum CMMProcessStep
{
	Start = 0,
	ChooseSettingsFinished = 1,
	ImageDetectedFinished = 2,
	PreOridentFinished = 3,
	BalanceFinished = 4,
	UnCodeMatcheFinished = 5,
	UseScaleBarsFinished = 6,
	EvaluationStep = 7,
};
static QString ReturnUnitStr(UnitType unit_type)
{
	switch (unit_type)
	{
	case Nanometer:
		return QObject::tr("nm");
		break;
	case Micrometer:
		return QApplication::translate("Micrometer", "\302\265m", 0);
		break;
	case Millimeter:
		return QObject::tr("mm");
		break;
	case Meter:
		return QObject::tr("m");
		break;
	case Kilometer:
		return QObject::tr("km");
		break;
	case Milligram:
		return QObject::tr("mg");
		break;
	case Gram:
		return QObject::tr("g");
		break;
	case Kilogram:
		return QObject::tr("kg");
		break;
	case Millisecond:
		return QObject::tr("ms");
		break;
	case Second:
		return QObject::tr("s");
		break;
	case Minute:
		return QObject::tr("min");
		break;
	case Hour:
		return QObject::tr("hour");
		break;
	case TempertureC:
		return QApplication::translate("TempertureC", "\342\204\203", 0);
		break;
	case TempertureF:
		return QApplication::translate("TempertureF", "\302\260F", 0);
		break;
	case TempertureK:
		return QObject::tr("K");
		break;
	case Nanoampere:
		return QObject::tr("nA");
		break;
	case Microampere:
		return QApplication::translate("Micrometer", "\302\265A", 0);
		break;
	case Milliampere:
		return QObject::tr("mA");
		break;
	case Ampere:
		return QObject::tr("A");
		break;
	case Kiloampere:
		return QObject::tr("kA");
		break;
	case Mole:
		return QObject::tr("mol");
		break;
	case Candela:
		return QObject::tr("cd");
		break;
	case Microstrain:
		return QApplication::translate("Microstrain", "\316\274\316\265", 0);
		break;
	case NoUnit:
		return QObject::tr("");
		break;
	case CUMM:
		return QApplication::translate("FittingElementDialog", "mm\302\263", 0);
		break;
	case CUCM:
		return QApplication::translate("FittingElementDialog", "cm\302\263", 0);
	case CUDM:
		return QApplication::translate("FittingElementDialog", "dm\302\263", 0);
	case CUM:
		return QApplication::translate("FittingElementDialog", "m\302\263", 0);
	}

	return QObject::tr("");
}

enum ImagePointType
{
	ImageReBuildCodePoint = 0,
	ImageReBuildUnCodePoint = 1,
	ImageUnReBuildCodePoint = 2,
	ImageUnReBuildUnCodePoint = 3
};
enum ObjectPointType
{
	ObjectReBuildCodePoint = 0,
	ObjectReBuildUnCodePoint = 1,
	ObjectUnReBuildCodePoint = 2,
	ObjectUnReBuildUnCodePoint = 3
};
enum Show3DDataType
{
	OnlyMorphology = 0,
	XData = 1,
	YData = 2,
	ZData = 3,
};
enum FitElementType
{
	FitEllipsoid = 0,
	FitPlane = 1,
};
#endif