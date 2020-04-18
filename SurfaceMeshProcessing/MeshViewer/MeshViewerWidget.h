#pragma once
#include <QString>
#include "QGLViewerWidget.h"
#include "MeshDefinition.h"

class MeshViewerWidget : public QGLViewerWidget
{
	Q_OBJECT
public:
	MeshViewerWidget(QWidget* parent = 0);
	virtual ~MeshViewerWidget(void);
	bool LoadMesh(const std::string & filename);
	void Clear(void);
	void UpdateMesh(void);
	bool SaveMesh(const std::string & filename);
	bool ScreenShot(void);
	void SetDrawBoundingBox(bool b);
	void SetDrawBoundary(bool b);
	//my
	void SetDrawBoundingSphere(bool b);
	void SetDrawVertexVector(bool b);
	//end my
	void EnableLighting(bool b);
	void EnableDoubleSide(bool b);
	void ResetView(void);
	void ViewCenter(void);
	void CopyRotation(void);
	void LoadRotation(void);
	//my function
	void PrintMeshAveArea(void);
	void PrintMeshAveCellCircum(void);
signals:
	void LoadMeshOKSignal(bool, QString);
public slots:
	void PrintMeshInfo(void);
protected:
	virtual void DrawScene(void) override;
	void DrawSceneMesh(void);

private:
	void DrawPoints(void) const;
	void DrawWireframe(void) const;
	void DrawHiddenLines(void) const;
	void DrawFlatLines(void) const;
	void DrawFlat(void) const;
	void DrawSmooth(void) const;
	void DrawBoundingBox(void) const;
	//my
	void DrawBoundingSphere(void) const;
	void DrawVertexNormalVector(void) const;
	void DrawSimpleBoundingSphere(void) const;
	void DrawMyPoint(double x, double y, double z, float r, float g, float b) const;
	void DrawMyLine(double x0, double y0, double z0, double x1, double y1, double z1) const;
	void DrawMyLine(double x0, double y0, double z0, double x1, double y1, double z1, float r, float g, float b) const;
	//end my
	void DrawBoundary(void) const;
protected:
	Mesh mesh;
	QString strMeshFileName;
	QString strMeshBaseName;
	QString strMeshPath;
	Mesh::Point ptMin;
	Mesh::Point ptMax;
	bool isEnableLighting;
	bool isTwoSideLighting;
	bool isDrawBoundingBox;
	bool isDrawBoundary;
	//my
	bool isDrawBoundingSphere;
	bool isDrawVertexVector;
	//end my
};
