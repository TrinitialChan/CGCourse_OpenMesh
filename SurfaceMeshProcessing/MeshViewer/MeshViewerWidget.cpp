#include <QtCore>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include "MeshViewerWidget.h"

#include <OpenMesh\Core\Mesh\TriMesh_ArrayKernelT.hh>

typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;

MeshViewerWidget::MeshViewerWidget(QWidget* parent)
	: QGLViewerWidget(parent),
	ptMin(0.0),
	ptMax(0.0),
	isEnableLighting(true),
	isTwoSideLighting(false),
	isDrawBoundingBox(false),
	isDrawBoundary(false),
	//my
	isDrawBoundingSphere(false),
	isDrawVertexVector(false)
	//end my
{
}

MeshViewerWidget::~MeshViewerWidget(void)
{

}

bool MeshViewerWidget::LoadMesh(const std::string & filename)
{
	Clear();
	bool read_OK = MeshTools::ReadMesh(mesh, filename);
	std::cout << "Load mesh from file " << filename << std::endl;
	if (read_OK)
	{
		strMeshFileName = QString::fromStdString(filename);
		QFileInfo fi(strMeshFileName);
		strMeshPath = fi.path();
		strMeshBaseName = fi.baseName();
		UpdateMesh();
		update();
		return true;
	}
	return false;
}

void MeshViewerWidget::Clear(void)
{
	mesh.clear();
}

void MeshViewerWidget::UpdateMesh(void)
{
	mesh.update_normals();
	if (mesh.vertices_empty())
	{
		std::cerr << "ERROR: UpdateMesh() No vertices!" << std::endl;
		return;
	}
	ptMin[0] = ptMin[1] = ptMin[2] = DBL_MAX;
	ptMax[0] = ptMax[1] = ptMax[2] = -DBL_MAX;
	for (const auto& vh : mesh.vertices())
	{
		ptMin.minimize(mesh.point(vh));
		ptMax.maximize(mesh.point(vh));
	}

	double avelen = 0.0;
	double maxlen = 0.0;
	double minlen = DBL_MAX;
	for (const auto& eh : mesh.edges())
	{
		double len = mesh.calc_edge_length(eh);
		maxlen = len > maxlen ? len : maxlen;
		minlen = len < minlen ? len : minlen;
		avelen += len;
	}

	SetScenePosition((ptMin + ptMax)*0.5, (ptMin - ptMax).norm()*0.5);
	std::cout << "Information of the input mesh:" << std::endl;
	std::cout << "  [V, E, F] = [" << mesh.n_vertices() << ", " << mesh.n_edges() << ", " << mesh.n_faces() << "]\n";
	std::cout << "  BoundingBox:\n";
	std::cout << "  X: [" << ptMin[0] << ", " << ptMax[0] << "]\n";
	std::cout << "  Y: [" << ptMin[1] << ", " << ptMax[1] << "]\n";
	std::cout << "  Z: [" << ptMin[2] << ", " << ptMax[2] << "]\n";
	std::cout << "  Diag length of BBox: " << (ptMax - ptMin).norm() << std::endl;
	std::cout << "  Edge Length: [" << minlen << ", " << maxlen << "]; AVG: " << avelen / mesh.n_edges() << std::endl;
}

bool MeshViewerWidget::SaveMesh(const std::string & filename)
{
	return MeshTools::WriteMesh(mesh, filename, DBL_DECIMAL_DIG);
}

bool MeshViewerWidget::ScreenShot()
{
	update();
	QString filename = strMeshPath + "/" + QDateTime::currentDateTime().toString("yyyyMMddHHmmsszzz") + QString(".png");
	QImage image = grabFramebuffer();
	image.save(filename);
	std::cout << "Save screen shot to " << filename.toStdString() << std::endl;
	return true;
}

void MeshViewerWidget::SetDrawBoundingBox(bool b)
{
	isDrawBoundingBox = b;
	update();
}
void MeshViewerWidget::SetDrawBoundary(bool b)
{
	isDrawBoundary = b;
	update();
}

//my
void MeshViewerWidget::SetDrawBoundingSphere(bool b)
{
	isDrawBoundingSphere = b;
	update();
}
void MeshViewerWidget::SetDrawVertexVector(bool b)
{
	isDrawVertexVector = b;
	update();
}
//end my

void MeshViewerWidget::EnableLighting(bool b)
{
	isEnableLighting = b;
	update();
}
void MeshViewerWidget::EnableDoubleSide(bool b)
{
	isTwoSideLighting = b;
	update();
}

void MeshViewerWidget::ResetView(void)
{
	ResetModelviewMatrix();
	ViewCenter();
	update();
}

void MeshViewerWidget::ViewCenter(void)
{
	if (!mesh.vertices_empty())
	{
		UpdateMesh();
	}
	update();
}

void MeshViewerWidget::CopyRotation(void)
{
	CopyModelViewMatrix();
}

void MeshViewerWidget::LoadRotation(void)
{
	LoadCopyModelViewMatrix();
	update();
}


//begin my
void MeshViewerWidget::PrintMeshAveArea(void)
{
	std::cout << "平均面积是\t"<<MeshTools::AveArea(mesh)<< std::endl;
	//std::cout << "shit";
}

void MeshViewerWidget::PrintMeshAveCellCircum(void)
{
	std::cout << "网格单元平均周长是\t"<<MeshTools::AveCellCircumlength(mesh) << std::endl;
	//std::cout << "shit";
}

//end my

void MeshViewerWidget::PrintMeshInfo(void)
{
	std::cout << "Mesh Info:\n";
	std::cout << "  [V, E, F] = [" << mesh.n_vertices() << ", " << mesh.n_edges() << ", " << mesh.n_faces() << "]\n";
	std::cout << "  BoundingBox:\n";
	std::cout << "  X: [" << ptMin[0] << ", " << ptMax[0] << "]\n";
	std::cout << "  Y: [" << ptMin[1] << ", " << ptMax[1] << "]\n";
	std::cout << "  Z: [" << ptMin[2] << ", " << ptMax[2] << "]\n";
	std::cout << "  Diag length of BBox: " << (ptMax - ptMin).norm() << std::endl;
}

void MeshViewerWidget::DrawScene(void)
{
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd(&projectionmatrix[0]);
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixd(&modelviewmatrix[0]);
	//DrawAxis();
	if (isDrawBoundingBox) DrawBoundingBox();
	//my
	if (isDrawBoundingSphere) DrawBoundingSphere();
	if (isDrawVertexVector) DrawVertexNormalVector();
	//end my

	if (isDrawBoundary) DrawBoundary();
	if (isEnableLighting) glEnable(GL_LIGHTING);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, isTwoSideLighting);
	DrawSceneMesh();
	if (isEnableLighting) glDisable(GL_LIGHTING);
}

void MeshViewerWidget::DrawSceneMesh(void)
{
	if (mesh.n_vertices() == 0) { return; }
	SetMaterial();
	switch (drawmode)
	{
	case POINTS:
		DrawPoints();
		break;
	case WIREFRAME:
		DrawWireframe();
		break;
	case HIDDENLINES:
		DrawHiddenLines();
		break;
	case FLATLINES:
		DrawFlatLines();
		break;
	case FLAT:
		glColor3d(0.8, 0.8, 0.8);
		DrawFlat();
		break;
	case SMOOTH:
		DrawSmooth();
		break;
	default:
		break;
	}
}

void MeshViewerWidget::DrawPoints(void) const
{
	glColor3d(1.0, 0.5, 0.5);
	glPointSize(5);
	glBegin(GL_POINTS);
	for (const auto& vh : mesh.vertices())
	{
		glNormal3dv(mesh.normal(vh).data());
		glVertex3dv(mesh.point(vh).data());
	}
	glEnd();
}

void MeshViewerWidget::DrawWireframe(void) const
{
	glColor3d(0.2, 0.2, 0.2);
	glBegin(GL_LINES);
	for (const auto& eh : mesh.edges())
	{
		auto heh = mesh.halfedge_handle(eh, 0);
		auto vh0 = mesh.from_vertex_handle(heh);
		auto vh1 = mesh.to_vertex_handle(heh);
		glNormal3dv(mesh.normal(vh0).data());
		glVertex3dv(mesh.point(vh0).data());
		glNormal3dv(mesh.normal(vh1).data());
		glVertex3dv(mesh.point(vh1).data());
	}
	glEnd();
}

void MeshViewerWidget::DrawHiddenLines() const
{
	glLineWidth(1.0);
	float backcolor[4];
	glGetFloatv(GL_COLOR_CLEAR_VALUE, backcolor);
	glColor4fv(backcolor);
	glDepthRange(0.01, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	if (glIsEnabled(GL_LIGHTING))
	{
		glDisable(GL_LIGHTING);
		DrawFlat();
		glEnable(GL_LIGHTING);
	}
	else
	{
		DrawFlat();
	}
	glDepthRange(0.0, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glColor3d(.3, .3, .3);
	DrawFlat();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void MeshViewerWidget::DrawFlatLines(void) const
{
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.5f, 2.0f);
	glShadeModel(GL_FLAT);
	//glColor3d(0.8, 0.8, 0.8);
	glColor3d(1.0, 1.0, 1.0);
	DrawFlat();
	glDisable(GL_POLYGON_OFFSET_FILL);
	if (glIsEnabled(GL_LIGHTING))
	{
		glDisable(GL_LIGHTING);
		DrawWireframe();
		glEnable(GL_LIGHTING);
	}
	else
	{
		DrawWireframe();
	}
}

void MeshViewerWidget::DrawFlat(void) const
{
	glBegin(GL_TRIANGLES);
	for (const auto& fh : mesh.faces())
	{
		glNormal3dv(mesh.normal(fh).data());
		for (const auto& fvh : mesh.fv_range(fh))
		{
			glVertex3dv(mesh.point(fvh).data());
		}
	}
	glEnd();
}

void MeshViewerWidget::DrawSmooth(void) const
{
	glColor3d(0.8, 0.8, 0.8);
	glShadeModel(GL_SMOOTH);
	glLoadName(static_cast<GLuint>(mesh.n_vertices()));
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_DOUBLE, 0, mesh.points());
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_DOUBLE, 0, mesh.vertex_normals());
	for (const auto& fh : mesh.faces())
	{
		glBegin(GL_POLYGON);
		for (const auto& fvh : mesh.fv_range(fh))
		{
			glArrayElement(fvh.idx());
		}
		glEnd();
	}
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
}

void MeshViewerWidget::DrawBoundingBox(void) const
{
	float linewidth;
	glGetFloatv(GL_LINE_WIDTH, &linewidth);
	glLineWidth(2.0f);
	glColor3d(.3, .7, .3);
	glBegin(GL_LINES);
	for (const auto& i : { 0, 1 })
	{
		for (const auto& j : { 0, 1 })
		{
			for (const auto& k : { 0, 1 })
			{
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(~i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], ~j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], ~k ? ptMin[2] : ptMax[2]);
			}
		}
	}
	glEnd();
	glLineWidth(linewidth);
	//DrawBoundingSphere();
	//DrawVertexNormalVector();
}
//my
void MeshViewerWidget::DrawSimpleBoundingSphere(void) const
{
	Mesh::VertexIter vlt, vBegin, vEnd;
	vBegin = mesh.vertices_begin();
	vEnd = mesh.vertices_end();
	double temp = 0;
	double maxd = 0;
	MeshTraits::Point maxpoint;
	MeshTraits::Point point;
	for (vlt = vBegin; vlt != vEnd; ++vlt) {
		point = mesh.point(*vlt);
		temp = sqrt(pow(point.data()[0], 2) + pow(point.data()[1], 2) + pow(point.data()[2], 2));
		if (temp > maxd) {
			maxd = temp;
			maxpoint = point;
		}
	}
	std::cout << "max_distance:" << maxd << "\n" << maxpoint.data()[0] << "\t" << maxpoint.data()[1] << "\t" << maxpoint.data()[2] << "\t" << std::endl;

	//绘制极大点
	glPointSize(10.0);
	glBegin(GL_POINTS);
	glColor3d(.3, .0, .0);
	glVertex3d(maxpoint.data()[0], maxpoint.data()[1], maxpoint.data()[2]);
	glEnd();

	glLineWidth(0.1f);
	glPointSize(1.0);
	glColor3d(.0, .0, .3);
	//glBegin(GL_LINE_STRIP);
	//glVertex3d(maxd, 0, 0);
	//glVertex3d(0, maxd, 0);
	//glVertex3d(0, 0, maxd);
	//double j = 0;

	//使用二次函数，采样效果差
	/*
	for (double j = -maxd; j <= maxd; j += 0.04) {
		glBegin(GL_LINE_STRIP);
		for (double i = -maxd; i <= maxd; i += 0.01)
		{
			glVertex3d(i, sqrt(pow(maxd, 2) - pow(i, 2)), j);
			glVertex3d(i, -sqrt(pow(maxd, 2) - pow(i, 2)), j);
		}
		glEnd();

		glBegin(GL_LINE_STRIP);
		for (double i = -maxd; i <= maxd; i += 0.01)
		{
			glVertex3d(j, sqrt(pow(maxd, 2) - pow(i, 2)), i);
			glVertex3d(j, -sqrt(pow(maxd, 2) - pow(i, 2)), i);
		}
		glEnd();

		glBegin(GL_LINE_STRIP);
		for (double i = -maxd; i <= maxd; i += 0.01)
		{
			glVertex3d(sqrt(pow(maxd, 2) - pow(i, 2)), j, i);
			glVertex3d(-sqrt(pow(maxd, 2) - pow(i, 2)), j, i);
		}
		glEnd();
	}
	*/
	double radius = maxd;
	//使用三角函数
	/*
	for (double j = -maxd; j <= maxd; j += 0.04) {
		glBegin(GL_LINE_STRIP);
		for (double i = 0; i <= 3.1415*2; i += 0.01)
		{
			glVertex3d(radius*cos(i), radius*sin(i), j);
		}
		glEnd();

		glBegin(GL_LINE_STRIP);
		for (double i = 0; i <= 3.1415 * 2; i += 0.01)
		{
			glVertex3d(j, radius*cos(i), radius*sin(i));
		}
		glEnd();

		glBegin(GL_LINE_STRIP);
		for (double i = 0; i <= 3.1415 * 2; i += 0.01)
		{
			glVertex3d(radius*sin(i), j, radius*cos(i));
		}
		glEnd();
	}
	*/
	//尝试绘制线框球体
	/*
	for (double j = -maxd; j <= maxd; j += 0.04) {
		radius = sqrt(maxd*maxd-j*j);

		glBegin(GL_LINE_STRIP);
		for (double i = 0; i <= 3.1415 * 2; i += 0.01)
		{
			glVertex3d(radius*cos(i), radius*sin(i), j);
		}
		glEnd();

		glBegin(GL_LINE_STRIP);
		for (double i = 0; i <= 3.1415 * 2; i += 0.01)
		{
			glVertex3d(j, radius*cos(i), radius*sin(i));
		}
		glEnd();

		glBegin(GL_LINE_STRIP);
		for (double i = 0; i <= 3.1415 * 2; i += 0.01)
		{
			glVertex3d(radius*sin(i), j, radius*cos(i));
		}
		glEnd();
	}*/
	//尝试绘制经纬线型
	for (double j = -maxd; j <= maxd; j += 0.04) {
		radius = sqrt(maxd*maxd - j * j);

		glBegin(GL_LINE_STRIP);
		for (double i = 0; i <= 3.1415 * 2; i += 0.01)
		{
			glVertex3d(radius*sin(i), j, radius*cos(i));
		}
		glEnd();
	}

	for (double j = 0; j < 3.1415; j += 0.03) {
		radius = maxd;
		glBegin(GL_LINE_STRIP);
		for (double i = 0; i <= 3.1415 * 2; i += 0.01)
		{
			glVertex3d(radius*cos(j)*cos(i), radius*sin(i), radius*sin(j)*cos(i));
		}
		glEnd();
	}
	
}
void MeshViewerWidget::DrawBoundingSphere(void) const
{
	Mesh::VertexIter vlt, vBegin, vEnd;
	vBegin = mesh.vertices_begin();
	vEnd = mesh.vertices_end();

	
	
	//初始化第一个球
	
	MeshTraits::Point maxpoint;
	MeshTraits::Point maxX,minX,maxY,minY,maxZ,minZ;
	double maxX_x = -1, maxY_y = -1, maxZ_z = -1;
	double minX_x = 1, minY_y = 1, minZ_z = 1;
	for (vlt = vBegin; vlt != vEnd; ++vlt) {
		MeshTraits::Point point = mesh.point(*vlt);

		if (point.data()[0] > maxX_x) { maxX_x = point.data()[0]; maxX = point; }
		if (point.data()[0] < minX_x) { minX_x = point.data()[0]; minX = point; }

		if (point.data()[1] > maxY_y) { maxY_y = point.data()[1]; maxY = point; }
		if (point.data()[1] < minY_y) { minY_y = point.data()[1]; minY = point; }

		if (point.data()[2] > maxZ_z) { maxZ_z = point.data()[2]; maxZ = point; }
		if (point.data()[2] < minZ_z) { minZ_z = point.data()[2]; minZ = point; }

	}

	//绘制极大点
	DrawMyPoint(maxX.data()[0], maxX.data()[1], maxX.data()[2], 0.8, 0, 0);
	DrawMyPoint(minX.data()[0], minX.data()[1], minX.data()[2], 0.8, 0, 0);

	DrawMyPoint(maxY.data()[0], maxY.data()[1], maxY.data()[2], 0.8, 0, 0);
	DrawMyPoint(minY.data()[0], minY.data()[1], minY.data()[2], 0.8, 0, 0);

	DrawMyPoint(maxZ.data()[0], maxZ.data()[1], maxZ.data()[2], 0.8, 0, 0);
	DrawMyPoint(minZ.data()[0], minZ.data()[1], minZ.data()[2], 0.8, 0, 0);

	double x=0, y=0, z=0;
	x = sqrt(pow((maxX.data()[0] - minX.data()[0]), 2) + pow((maxX.data()[1] - minX.data()[1]), 2) + pow((maxX.data()[2] - minX.data()[2]), 2));
	y = sqrt(pow((maxY.data()[0] - minY.data()[0]), 2) + pow((maxY.data()[1] - minY.data()[1]), 2) + pow((maxY.data()[2] - minY.data()[2]), 2));
	z = sqrt(pow((maxZ.data()[0] - minZ.data()[0]), 2) + pow((maxZ.data()[1] - minZ.data()[1]), 2) + pow((maxZ.data()[2] - minZ.data()[2]), 2));

	float dia = 0;
	MeshTraits::Point max = maxX, min = minX;
	dia = x;
	if (z > x && z > y)
	{
		max = maxZ;
		min = minZ;
		dia = z;
	}
	else if (y > x && y > z)
	{
		max = maxY;
		min = minY;
		dia = y;
	}

	//初始球心
	double tempx, tempy, tempz;
	double tempRadius;
	tempx = 0.5*(max.data()[0] + min.data()[0]);
	tempy = 0.5*(max.data()[1] + min.data()[1]);
	tempz = 0.5*(max.data()[2] + min.data()[2]);
	
	//初始球半径
	tempRadius = dia/2;

	//绘制初始球心
	std::cout << "初始球心是:\n" << tempx << "\t" << tempy << "\t" << tempz << std::endl;
	DrawMyPoint(tempx, tempy, tempz, 1, 0.3, 0.5);
	std::cout << "初始半径是:\t" << tempRadius << std::endl;

	//开始修正
	for (vlt = vBegin; vlt != vEnd; ++vlt) {
		MeshTraits::Point point = mesh.point(*vlt);
		double x, y, z;
		x = point.data()[0]; y = point.data()[1]; z = point.data()[2];
		double dist = sqrt(pow((tempx - x), 2) + pow((tempy - y), 2) + pow((tempz - z), 2));

		if (dist > tempRadius) {
			DrawMyPoint(x, y, z,0.8,0,0.8);
			double bias = (dist - tempRadius)/2;
			tempx += (x - tempx)*bias / dist;
			tempy += (y - tempy)*bias / dist;
			tempz += (z - tempz)*bias / dist;
			tempRadius = (dist + tempRadius) / 2;
		}

	}
	//修正完成，标注球心
	DrawMyPoint(tempx, tempy, tempz,0.8,0.8,0);
	std::cout <<"修正球心是:\n"<< tempx << "\t"<<tempy << "\t" << tempz << std::endl;
	std::cout << "修正半径是:\t" << tempRadius << std::endl;

	//定义半径
	double radius;
	
	//定义球心
	double center_x, center_y, center_z=0;
	
	//传参
	center_x = tempx;
	center_y = tempy;
	center_z = tempz;
	radius = tempRadius;

	//绘制经纬线型
	glLineWidth(0.1f);
	glPointSize(1.0);
	glColor3d(.0, .0, .3);
	for (double j = -radius; j <= radius; j += 0.04) {
		double layer_radius = sqrt(radius*radius - j * j);

		glBegin(GL_LINE_STRIP);
		for (double i = 0; i <= 3.1415 * 2; i += 0.01)
		{
			glVertex3d(	layer_radius*sin(i)+ center_x ,
						j + center_y, 
						layer_radius*cos(i) + center_z);
		}
		glEnd();
	}

	for (double j = 0; j < 3.1415; j += 0.03) {
		glBegin(GL_LINE_STRIP);
		for (double i = 0; i <= 3.1415 * 2; i += 0.01)
		{
			glVertex3d(	radius*cos(j)*cos(i) + center_x, 
						radius*sin(i) + center_y, 
						radius*sin(j)*cos(i) + center_z);
		}
		glEnd();
	}

}
void MeshViewerWidget::DrawVertexNormalVector(void) const 
{

	Mesh::VertexIter vlt, vBegin, vEnd;
	vBegin = mesh.vertices_begin();
	vEnd = mesh.vertices_end();

	for (vlt = vBegin; vlt != vEnd; ++vlt) {
		MeshTraits::Point point = mesh.point(*vlt);
		//auto vertex = mesh.vertex(*vlt);
		OpenMesh::VertexHandle vertex = *vlt;
		
		//for (auto it = mesh.cvf_begin(vertex); it != mesh.cvf_end(vertex); ++it) {}

		Mesh::VertexFaceIter vflt, vfBegin, vfEnd;
		vfBegin = mesh.cvf_begin(*vlt);
		vfEnd = mesh.cvf_end(*vlt);
		double sum = 0;
		double dotx=0, doty=0, dotz=0;
		for (vflt= vfBegin; vflt != vfEnd; ++vflt) {

			//auto face = vflt.handle();
			//MyMesh::Face face = mesh.face(*vflt);
			auto normal = mesh.normal(*vflt);
			double x = normal.data()[0];
			double y = normal.data()[1];
			double z = normal.data()[2];
			
			dotx += x;
			doty += y;
			dotz += z;
			sum++;
		}
		double myMod = sqrt(dotx * dotx + doty * doty + dotz * dotz);
		myMod;
		dotx /= myMod;
		doty /= myMod;
		dotz /= myMod;
		double ratio = 0.01;
		dotx *= ratio;
		doty *= ratio;
		dotz *= ratio;
		DrawMyLine(point.data()[0], point.data()[1], point.data()[2], point.data()[0] + dotx, point.data()[1] + doty, point.data()[2] + dotz);
		DrawMyLine(point.data()[0], point.data()[1], point.data()[2], point.data()[0] - dotx, point.data()[1] - doty, point.data()[2] - dotz);
		//std::cout << sum << std::endl;

	}

}
void MeshViewerWidget::DrawMyPoint(double x,double y, double z,float r,float g,float b) const {
	glPointSize(10.0);
	glBegin(GL_POINTS);
	glColor3d(r, g, b);
	glVertex3d(x,y,z);
	glEnd();
}

void MeshViewerWidget::DrawMyLine(double x0, double y0, double z0, double x1, double y1, double z1) const 
{
	DrawMyLine(x0, y0, z0, x1, y1, z1, 0.3, 0, 0);
}

void MeshViewerWidget::DrawMyLine(double x0, double y0, double z0, double x1, double y1, double z1, float r, float g, float b) const
{
	glLineWidth(0.1f);
	glBegin(GL_LINE_STRIP);
	glColor3d(r, g, b);
	glVertex3d(x0, y0, z0);
	glVertex3d(x1, y1, z1);
	glEnd();
}
//end my
void MeshViewerWidget::DrawBoundary(void) const
{
	float linewidth;
	glGetFloatv(GL_LINE_WIDTH, &linewidth);
	glLineWidth(2.0f);
	glColor3d(0.1, 0.1, 0.1);
	glBegin(GL_LINES);
	for (const auto& eh : mesh.edges())
	{
		if (mesh.is_boundary(eh))
		{
			auto heh = mesh.halfedge_handle(eh, 0);
			auto vh0 = mesh.from_vertex_handle(heh);
			auto vh1 = mesh.to_vertex_handle(heh);
			glNormal3dv(mesh.normal(vh0).data());
			glVertex3dv(mesh.point(vh0).data());
			glNormal3dv(mesh.normal(vh1).data());
			glVertex3dv(mesh.point(vh1).data());
		}
	}
	glEnd();
	glLineWidth(linewidth);
}

