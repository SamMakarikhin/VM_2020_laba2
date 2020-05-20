#pragma once
#include<math.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
# define PI           3.14159265358979323846

//using namespace std;
namespace Laba2 {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;
	double a = 1, l = 7, h = 0.1, tau = 1, T = 3;
	double B1, B2, B3, B4, B5;
	double **y;

	double phi(int a1, int a2, int a3, double x) {
		return a1 * 1 + a2*cos(x*PI / l) + a3*cos(2 * PI*x / l);
	}
	double psi(int a1, int a2, int a3, double x) {
		return a1 * 1 + a2*cos(x*PI / l) + a3*cos(2 * PI*x / l);
	}
	double bx(double x) {
		return B1 * 1 + B2*cos(x*PI / l) + B3*sin(x*PI / l) + B4*cos(2 * PI*x / l) + B5*cos(2 * PI*x / l);
	}

	double Simpson(double(*f)(double), double l){
		double h = l / 1000;
		double S4 = f(0 + h), S2 = 0.;
		for (int i = 3; i < 1000; i += 2) {
			S4 += f(0 + i*h);
			S2 += f(0 + (i - 1)*h);
		}
		return (h / 3) * (f(0) + 4 * S4 + 2 * S2 + f(l));
	}

	double** GenerateMatrix(double a, double tau, double h, int sizeN){
		double** A = new double*[sizeN];
		for (int i = 0; i < sizeN; i++) {
			A[i] = new double[sizeN];
		}

		for (int i = 0; i < sizeN; i++){
			for (int j = 0; j < sizeN; j++) {
				if (i == j) A[i][j] = (1 / (tau*tau)) - (2 / (h*h));
				else if (j == i + 1) A[i][j] = (-1)*a / (h*h);
				else if (j == i - 1) A[i][j] = a / (h*h);
				else A[i][j] = 0;
			}
		}
		return A;
	}


	double* GenerateVector(double** y, double(*f)(double),double tau, double h, int sizeN, int j){
		double* vec = new double[sizeN];
		for (int i = 0; i < sizeN; i++){
			if (i == 0 || i == sizeN - 1) vec[i] = 0;
			else
			vec[i] = 2 * y[i][j-1] / (tau*tau) - y[i][j-2] / (tau*tau) + y[i][j]* ( bx(i*h) + y[i][0] * Simpson(&bx,l));
		}
		return vec;
			
	}

/*	double* diagonal(double** A, double* b, int N){
		double y;
		double *a, *bb, *res;
		a = new double[N];
		bb = new double[N];
		res = new double[N];
		y = A[0][0];
		a[0] = -A[0][1] / y;
		bb[0] = b[0] / y;
		for (int i = 1; i < N-1; i++) {
			y = A[i][i] + A[i][i - 1] * a[i - 1];
			a[i] = -A[i][i + 1] / y;
			bb[i] = (b[i] - A[i][i - 1] * b[i - 1]) / y;
		}
	}*/
    def sweep_method(a, b, c, func, count)
    {
      A = []
      B = []
      res = [0] * count
      A.append(-c[0] / b[0])
      B.append(func[0] / b[0])
      for i in range(1, count)
      {
        A.append(-c[i] / (a[i] * A[i - 1] + b[i]))
        B.append((func[i] - a[i] * B[i - 1]) / (a[i] * A[i - 1] + b[i]))
        res[count - 1] = B[count - 1]
      }
      for i in range(count - 2, -1, -1)
      {
        res[i] = (A[i] * res[i + 1] + B[i])      }      return res    }
	double * gauss(double **a, double *y, int n)
	{
		double *x, max;
		int k, index;
		const double eps = 0.00001;  // точность
		x = new double[n];
		k = 0;
		while (k < n)
		{
			// Поиск строки с максимальным a[i][k]
			max = abs(a[k][k]);
			index = k;
			for (int i = k + 1; i < n; i++)
			{
				if (abs(a[i][k]) > max)
				{
					max = abs(a[i][k]);
					index = i;
				}
			}
			// Перестановка строк
			if (max < eps)
			{
				// нет ненулевых диагональных элементов
				std::cout << "Решение получить невозможно из-за нулевого столбца ";
				std::cout << index << " матрицы A" << std::endl;
				return 0;
			}
			for (int j = 0; j < n; j++)
			{
				double temp = a[k][j];
				a[k][j] = a[index][j];
				a[index][j] = temp;
			}
			double temp = y[k];
			y[k] = y[index];
			y[index] = temp;
			// Нормализация уравнений
			for (int i = k; i < n; i++)
			{
				double temp = a[i][k];
				if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
				for (int j = 0; j < n; j++)
					a[i][j] = a[i][j] / temp;
				y[i] = y[i] / temp;
				if (i == k)  continue; // уравнение не вычитать само из себя
				for (int j = 0; j < n; j++)
					a[i][j] = a[i][j] - a[k][j];
				y[i] = y[i] - y[k];
			}
			k++;
		}
		// обратная подстановка
		for (k = n - 1; k >= 0; k--)
		{
			x[k] = y[k];
			for (int i = 0; i < k; i++)
				y[i] = y[i] - a[i][k] * x[k];
		}
		return x;
	}

	/// <summary>
	/// Summary for MyForm
	/// </summary>
	public ref class MyForm : public System::Windows::Forms::Form
	{
	public:
		MyForm(void)
		{
			InitializeComponent();
			//
			//TODO: Add the constructor code here
			//
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~MyForm()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::Label^  label1;
	protected:
	private: System::Windows::Forms::Label^  label2;
	private: System::Windows::Forms::TextBox^  textBox1;
	private: System::Windows::Forms::TextBox^  textBox2;
	private: System::Windows::Forms::TextBox^  textBox3;
	private: System::Windows::Forms::TextBox^  textBox4;
	private: System::Windows::Forms::TextBox^  textBox5;
	private: System::Windows::Forms::TextBox^  textBox6;
	private: System::Windows::Forms::Label^  label3;
	private: System::Windows::Forms::Label^  label4;
	private: System::Windows::Forms::Label^  label5;
	private: System::Windows::Forms::Label^  label6;
	private: System::Windows::Forms::Label^  label7;
	private: System::Windows::Forms::Label^  label8;
	private: System::Windows::Forms::DataVisualization::Charting::Chart^  chart1;
	private: System::Windows::Forms::Button^  button1;
	private: System::Windows::Forms::Label^  label9;
	private: System::Windows::Forms::TextBox^  textBox7;
	private: System::Windows::Forms::TextBox^  textBox8;
	private: System::Windows::Forms::TextBox^  textBox9;
	private: System::Windows::Forms::TextBox^  textBox10;
	private: System::Windows::Forms::TextBox^  textBox11;
	private: System::Windows::Forms::Label^  label10;
	private: System::Windows::Forms::Label^  label11;
	private: System::Windows::Forms::Label^  label12;
	private: System::Windows::Forms::Label^  label13;
	private: System::Windows::Forms::Label^  label14;
	private: System::Windows::Forms::Label^  label15;
	private: System::Windows::Forms::TextBox^  textBox12;
	private: System::Windows::Forms::TextBox^  textBox13;
	private: System::Windows::Forms::TextBox^  textBox14;
	private: System::Windows::Forms::TextBox^  textBox15;
	private: System::Windows::Forms::TextBox^  textBox16;
	private: System::Windows::Forms::Label^  label16;
	private: System::Windows::Forms::Label^  label17;
	private: System::Windows::Forms::Label^  label18;
	private: System::Windows::Forms::Label^  label19;
	private: System::Windows::Forms::Label^  label20;
	private: System::Windows::Forms::Label^  label21;
	private: System::Windows::Forms::Label^  label22;
	private: System::Windows::Forms::Label^  label23;

	private: System::ComponentModel::IContainer^  components;




	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>


#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			System::Windows::Forms::DataVisualization::Charting::ChartArea^  chartArea5 = (gcnew System::Windows::Forms::DataVisualization::Charting::ChartArea());
			System::Windows::Forms::DataVisualization::Charting::Legend^  legend5 = (gcnew System::Windows::Forms::DataVisualization::Charting::Legend());
			System::Windows::Forms::DataVisualization::Charting::Series^  series9 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			System::Windows::Forms::DataVisualization::Charting::Series^  series10 = (gcnew System::Windows::Forms::DataVisualization::Charting::Series());
			System::Windows::Forms::DataVisualization::Charting::Title^  title5 = (gcnew System::Windows::Forms::DataVisualization::Charting::Title());
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->textBox1 = (gcnew System::Windows::Forms::TextBox());
			this->textBox2 = (gcnew System::Windows::Forms::TextBox());
			this->textBox3 = (gcnew System::Windows::Forms::TextBox());
			this->textBox4 = (gcnew System::Windows::Forms::TextBox());
			this->textBox5 = (gcnew System::Windows::Forms::TextBox());
			this->textBox6 = (gcnew System::Windows::Forms::TextBox());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->label4 = (gcnew System::Windows::Forms::Label());
			this->label5 = (gcnew System::Windows::Forms::Label());
			this->label6 = (gcnew System::Windows::Forms::Label());
			this->label7 = (gcnew System::Windows::Forms::Label());
			this->label8 = (gcnew System::Windows::Forms::Label());
			this->chart1 = (gcnew System::Windows::Forms::DataVisualization::Charting::Chart());
			this->button1 = (gcnew System::Windows::Forms::Button());
			this->label9 = (gcnew System::Windows::Forms::Label());
			this->textBox7 = (gcnew System::Windows::Forms::TextBox());
			this->textBox8 = (gcnew System::Windows::Forms::TextBox());
			this->textBox9 = (gcnew System::Windows::Forms::TextBox());
			this->textBox10 = (gcnew System::Windows::Forms::TextBox());
			this->textBox11 = (gcnew System::Windows::Forms::TextBox());
			this->label10 = (gcnew System::Windows::Forms::Label());
			this->label11 = (gcnew System::Windows::Forms::Label());
			this->label12 = (gcnew System::Windows::Forms::Label());
			this->label13 = (gcnew System::Windows::Forms::Label());
			this->label14 = (gcnew System::Windows::Forms::Label());
			this->label15 = (gcnew System::Windows::Forms::Label());
			this->textBox12 = (gcnew System::Windows::Forms::TextBox());
			this->textBox13 = (gcnew System::Windows::Forms::TextBox());
			this->textBox14 = (gcnew System::Windows::Forms::TextBox());
			this->textBox15 = (gcnew System::Windows::Forms::TextBox());
			this->textBox16 = (gcnew System::Windows::Forms::TextBox());
			this->label16 = (gcnew System::Windows::Forms::Label());
			this->label17 = (gcnew System::Windows::Forms::Label());
			this->label18 = (gcnew System::Windows::Forms::Label());
			this->label19 = (gcnew System::Windows::Forms::Label());
			this->label20 = (gcnew System::Windows::Forms::Label());
			this->label21 = (gcnew System::Windows::Forms::Label());
			this->label22 = (gcnew System::Windows::Forms::Label());
			this->label23 = (gcnew System::Windows::Forms::Label());
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->chart1))->BeginInit();
			this->SuspendLayout();
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label1->Location = System::Drawing::Point(8, 46);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(48, 19);
			this->label1->TabIndex = 0;
			this->label1->Text = L"φ(х) =";
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label2->Location = System::Drawing::Point(8, 90);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(49, 19);
			this->label2->TabIndex = 1;
			this->label2->Text = L"ψ(х) =";
			// 
			// textBox1
			// 
			this->textBox1->Location = System::Drawing::Point(62, 47);
			this->textBox1->Name = L"textBox1";
			this->textBox1->Size = System::Drawing::Size(43, 20);
			this->textBox1->TabIndex = 2;
			this->textBox1->Text = L"1";
			this->textBox1->TextChanged += gcnew System::EventHandler(this, &MyForm::textBox1_TextChanged);
			// 
			// textBox2
			// 
			this->textBox2->Location = System::Drawing::Point(127, 48);
			this->textBox2->Name = L"textBox2";
			this->textBox2->Size = System::Drawing::Size(43, 20);
			this->textBox2->TabIndex = 3;
			this->textBox2->Text = L"1";
			// 
			// textBox3
			// 
			this->textBox3->Location = System::Drawing::Point(236, 50);
			this->textBox3->Name = L"textBox3";
			this->textBox3->Size = System::Drawing::Size(43, 20);
			this->textBox3->TabIndex = 4;
			this->textBox3->Text = L"1";
			// 
			// textBox4
			// 
			this->textBox4->Location = System::Drawing::Point(62, 91);
			this->textBox4->Name = L"textBox4";
			this->textBox4->Size = System::Drawing::Size(43, 20);
			this->textBox4->TabIndex = 5;
			this->textBox4->Text = L"1";
			// 
			// textBox5
			// 
			this->textBox5->Location = System::Drawing::Point(127, 90);
			this->textBox5->Name = L"textBox5";
			this->textBox5->Size = System::Drawing::Size(43, 20);
			this->textBox5->TabIndex = 6;
			this->textBox5->Text = L"1";
			// 
			// textBox6
			// 
			this->textBox6->Location = System::Drawing::Point(236, 91);
			this->textBox6->Name = L"textBox6";
			this->textBox6->Size = System::Drawing::Size(43, 20);
			this->textBox6->TabIndex = 7;
			this->textBox6->Text = L"1";
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label3->Location = System::Drawing::Point(107, 47);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(18, 19);
			this->label3->TabIndex = 8;
			this->label3->Text = L"+";
			// 
			// label4
			// 
			this->label4->AutoSize = true;
			this->label4->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label4->Location = System::Drawing::Point(173, 49);
			this->label4->Name = L"label4";
			this->label4->Size = System::Drawing::Size(63, 19);
			this->label4->TabIndex = 9;
			this->label4->Text = L"*cos α +";
			// 
			// label5
			// 
			this->label5->AutoSize = true;
			this->label5->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label5->Location = System::Drawing::Point(279, 50);
			this->label5->Name = L"label5";
			this->label5->Size = System::Drawing::Size(50, 19);
			this->label5->TabIndex = 10;
			this->label5->Text = L"*cos β";
			// 
			// label6
			// 
			this->label6->AutoSize = true;
			this->label6->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label6->Location = System::Drawing::Point(107, 91);
			this->label6->Name = L"label6";
			this->label6->Size = System::Drawing::Size(18, 19);
			this->label6->TabIndex = 11;
			this->label6->Text = L"+";
			// 
			// label7
			// 
			this->label7->AutoSize = true;
			this->label7->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label7->Location = System::Drawing::Point(173, 91);
			this->label7->Name = L"label7";
			this->label7->Size = System::Drawing::Size(63, 19);
			this->label7->TabIndex = 12;
			this->label7->Text = L"*cos α +";
			// 
			// label8
			// 
			this->label8->AutoSize = true;
			this->label8->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label8->Location = System::Drawing::Point(279, 91);
			this->label8->Name = L"label8";
			this->label8->Size = System::Drawing::Size(50, 19);
			this->label8->TabIndex = 13;
			this->label8->Text = L"*cos β";
			// 
			// chart1
			// 
			chartArea5->Name = L"ChartArea1";
			this->chart1->ChartAreas->Add(chartArea5);
			legend5->Name = L"Legend1";
			this->chart1->Legends->Add(legend5);
			this->chart1->Location = System::Drawing::Point(12, 183);
			this->chart1->Name = L"chart1";
			this->chart1->Palette = System::Windows::Forms::DataVisualization::Charting::ChartColorPalette::Berry;
			series9->ChartArea = L"ChartArea1";
			series9->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Spline;
			series9->Color = System::Drawing::Color::Blue;
			series9->Legend = L"Legend1";
			series9->Name = L"Series1";
			series10->BorderWidth = 2;
			series10->ChartArea = L"ChartArea1";
			series10->ChartType = System::Windows::Forms::DataVisualization::Charting::SeriesChartType::Spline;
			series10->Color = System::Drawing::Color::Red;
			series10->Legend = L"Legend1";
			series10->Name = L"Series2";
			this->chart1->Series->Add(series9);
			this->chart1->Series->Add(series10);
			this->chart1->Size = System::Drawing::Size(899, 376);
			this->chart1->TabIndex = 14;
			this->chart1->Text = L"chart1";
			title5->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 12));
			title5->Name = L"Title1";
			title5->Text = L" ";
			this->chart1->Titles->Add(title5);
			this->chart1->Click += gcnew System::EventHandler(this, &MyForm::chart1_Click);
			// 
			// button1
			// 
			this->button1->Font = (gcnew System::Drawing::Font(L"Times New Roman", 15.75F, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->button1->Location = System::Drawing::Point(775, 12);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(136, 145);
			this->button1->TabIndex = 15;
			this->button1->Text = L"ИТОГ";
			this->button1->UseVisualStyleBackColor = true;
			this->button1->Click += gcnew System::EventHandler(this, &MyForm::button1_Click);
			// 
			// label9
			// 
			this->label9->AutoSize = true;
			this->label9->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label9->Location = System::Drawing::Point(8, 134);
			this->label9->Name = L"label9";
			this->label9->Size = System::Drawing::Size(47, 19);
			this->label9->TabIndex = 16;
			this->label9->Text = L"b(x) =";
			// 
			// textBox7
			// 
			this->textBox7->Location = System::Drawing::Point(62, 135);
			this->textBox7->Name = L"textBox7";
			this->textBox7->Size = System::Drawing::Size(43, 20);
			this->textBox7->TabIndex = 17;
			this->textBox7->Text = L"1";
			// 
			// textBox8
			// 
			this->textBox8->Location = System::Drawing::Point(127, 134);
			this->textBox8->Name = L"textBox8";
			this->textBox8->Size = System::Drawing::Size(43, 20);
			this->textBox8->TabIndex = 18;
			this->textBox8->Text = L"1";
			// 
			// textBox9
			// 
			this->textBox9->Location = System::Drawing::Point(236, 135);
			this->textBox9->Name = L"textBox9";
			this->textBox9->Size = System::Drawing::Size(43, 20);
			this->textBox9->TabIndex = 19;
			this->textBox9->Text = L"1";
			// 
			// textBox10
			// 
			this->textBox10->Location = System::Drawing::Point(338, 137);
			this->textBox10->Name = L"textBox10";
			this->textBox10->Size = System::Drawing::Size(43, 20);
			this->textBox10->TabIndex = 20;
			this->textBox10->Text = L"1";
			// 
			// textBox11
			// 
			this->textBox11->Location = System::Drawing::Point(445, 137);
			this->textBox11->Name = L"textBox11";
			this->textBox11->Size = System::Drawing::Size(43, 20);
			this->textBox11->TabIndex = 21;
			this->textBox11->Text = L"1";
			// 
			// label10
			// 
			this->label10->AutoSize = true;
			this->label10->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label10->Location = System::Drawing::Point(556, 75);
			this->label10->Name = L"label10";
			this->label10->Size = System::Drawing::Size(57, 19);
			this->label10->TabIndex = 22;
			this->label10->Text = L"a ^ 2 = ";
			// 
			// label11
			// 
			this->label11->AutoSize = true;
			this->label11->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label11->Location = System::Drawing::Point(107, 136);
			this->label11->Name = L"label11";
			this->label11->Size = System::Drawing::Size(18, 19);
			this->label11->TabIndex = 23;
			this->label11->Text = L"+";
			// 
			// label12
			// 
			this->label12->AutoSize = true;
			this->label12->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label12->Location = System::Drawing::Point(173, 135);
			this->label12->Name = L"label12";
			this->label12->Size = System::Drawing::Size(63, 19);
			this->label12->TabIndex = 24;
			this->label12->Text = L"*cos α +";
			// 
			// label13
			// 
			this->label13->AutoSize = true;
			this->label13->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label13->Location = System::Drawing::Point(279, 136);
			this->label13->Name = L"label13";
			this->label13->Size = System::Drawing::Size(58, 19);
			this->label13->TabIndex = 25;
			this->label13->Text = L"*sin α +";
			this->label13->Click += gcnew System::EventHandler(this, &MyForm::label13_Click);
			// 
			// label14
			// 
			this->label14->AutoSize = true;
			this->label14->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label14->Location = System::Drawing::Point(381, 136);
			this->label14->Name = L"label14";
			this->label14->Size = System::Drawing::Size(63, 19);
			this->label14->TabIndex = 26;
			this->label14->Text = L"*cos β +";
			// 
			// label15
			// 
			this->label15->AutoSize = true;
			this->label15->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label15->Location = System::Drawing::Point(489, 138);
			this->label15->Name = L"label15";
			this->label15->Size = System::Drawing::Size(45, 19);
			this->label15->TabIndex = 27;
			this->label15->Text = L"*sin β";
			// 
			// textBox12
			// 
			this->textBox12->Location = System::Drawing::Point(619, 76);
			this->textBox12->Name = L"textBox12";
			this->textBox12->Size = System::Drawing::Size(38, 20);
			this->textBox12->TabIndex = 28;
			this->textBox12->Text = L"1";
			// 
			// textBox13
			// 
			this->textBox13->Location = System::Drawing::Point(717, 109);
			this->textBox13->Name = L"textBox13";
			this->textBox13->Size = System::Drawing::Size(43, 20);
			this->textBox13->TabIndex = 29;
			// 
			// textBox14
			// 
			this->textBox14->Location = System::Drawing::Point(619, 105);
			this->textBox14->Name = L"textBox14";
			this->textBox14->Size = System::Drawing::Size(38, 20);
			this->textBox14->TabIndex = 30;
			this->textBox14->Text = L"7";
			// 
			// textBox15
			// 
			this->textBox15->Location = System::Drawing::Point(619, 137);
			this->textBox15->Name = L"textBox15";
			this->textBox15->Size = System::Drawing::Size(38, 20);
			this->textBox15->TabIndex = 31;
			this->textBox15->Text = L"0,1";
			// 
			// textBox16
			// 
			this->textBox16->Location = System::Drawing::Point(717, 137);
			this->textBox16->Name = L"textBox16";
			this->textBox16->Size = System::Drawing::Size(43, 20);
			this->textBox16->TabIndex = 32;
			// 
			// label16
			// 
			this->label16->AutoSize = true;
			this->label16->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label16->Location = System::Drawing::Point(557, 104);
			this->label16->Name = L"label16";
			this->label16->Size = System::Drawing::Size(29, 19);
			this->label16->TabIndex = 33;
			this->label16->Text = L"l = ";
			// 
			// label17
			// 
			this->label17->AutoSize = true;
			this->label17->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label17->Location = System::Drawing::Point(680, 110);
			this->label17->Name = L"label17";
			this->label17->Size = System::Drawing::Size(31, 19);
			this->label17->TabIndex = 34;
			this->label17->Text = L"T =";
			// 
			// label18
			// 
			this->label18->AutoSize = true;
			this->label18->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label18->Location = System::Drawing::Point(557, 138);
			this->label18->Name = L"label18";
			this->label18->Size = System::Drawing::Size(29, 19);
			this->label18->TabIndex = 35;
			this->label18->Text = L"h =";
			// 
			// label19
			// 
			this->label19->AutoSize = true;
			this->label19->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label19->Location = System::Drawing::Point(682, 136);
			this->label19->Name = L"label19";
			this->label19->Size = System::Drawing::Size(30, 19);
			this->label19->TabIndex = 36;
			this->label19->Text = L"t = ";
			// 
			// label20
			// 
			this->label20->AutoSize = true;
			this->label20->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Bold, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label20->Location = System::Drawing::Point(8, 10);
			this->label20->Name = L"label20";
			this->label20->Size = System::Drawing::Size(134, 19);
			this->label20->TabIndex = 37;
			this->label20->Text = L"Введите данные:";
			// 
			// label21
			// 
			this->label21->AutoSize = true;
			this->label21->Font = (gcnew System::Drawing::Font(L"Times New Roman", 15.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label21->Location = System::Drawing::Point(369, 93);
			this->label21->Name = L"label21";
			this->label21->Size = System::Drawing::Size(122, 23);
			this->label21->TabIndex = 38;
			this->label21->Text = L"β = 2*π*k / 7";
			// 
			// label22
			// 
			this->label22->AutoSize = true;
			this->label22->Font = (gcnew System::Drawing::Font(L"Times New Roman", 15.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label22->Location = System::Drawing::Point(369, 54);
			this->label22->Name = L"label22";
			this->label22->Size = System::Drawing::Size(101, 23);
			this->label22->TabIndex = 39;
			this->label22->Text = L"α = π*k / 7";
			// 
			// label23
			// 
			this->label23->AutoSize = true;
			this->label23->Font = (gcnew System::Drawing::Font(L"Microsoft Sans Serif", 15.75F, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->label23->Location = System::Drawing::Point(524, 31);
			this->label23->Name = L"label23";
			this->label23->Size = System::Drawing::Size(0, 25);
			this->label23->TabIndex = 40;
			// 
			// MyForm
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(923, 571);
			this->Controls->Add(this->label23);
			this->Controls->Add(this->label22);
			this->Controls->Add(this->label21);
			this->Controls->Add(this->label20);
			this->Controls->Add(this->textBox16);
			this->Controls->Add(this->label10);
			this->Controls->Add(this->label19);
			this->Controls->Add(this->textBox13);
			this->Controls->Add(this->label15);
			this->Controls->Add(this->label17);
			this->Controls->Add(this->textBox15);
			this->Controls->Add(this->label14);
			this->Controls->Add(this->textBox14);
			this->Controls->Add(this->label18);
			this->Controls->Add(this->label13);
			this->Controls->Add(this->label12);
			this->Controls->Add(this->textBox12);
			this->Controls->Add(this->label16);
			this->Controls->Add(this->label11);
			this->Controls->Add(this->textBox11);
			this->Controls->Add(this->textBox10);
			this->Controls->Add(this->textBox9);
			this->Controls->Add(this->textBox8);
			this->Controls->Add(this->textBox7);
			this->Controls->Add(this->label9);
			this->Controls->Add(this->button1);
			this->Controls->Add(this->chart1);
			this->Controls->Add(this->label8);
			this->Controls->Add(this->label7);
			this->Controls->Add(this->label6);
			this->Controls->Add(this->label5);
			this->Controls->Add(this->label4);
			this->Controls->Add(this->label3);
			this->Controls->Add(this->textBox6);
			this->Controls->Add(this->textBox5);
			this->Controls->Add(this->textBox4);
			this->Controls->Add(this->textBox3);
			this->Controls->Add(this->textBox2);
			this->Controls->Add(this->textBox1);
			this->Controls->Add(this->label2);
			this->Controls->Add(this->label1);
			this->Name = L"MyForm";
			this->Text = L"MyForm";
			this->Load += gcnew System::EventHandler(this, &MyForm::MyForm_Load);
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->chart1))->EndInit();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
	private: System::Void button1_Click(System::Object^  sender, System::EventArgs^  e) {
				 chart1->Series["Series1"]->Points->Clear();
				 chart1->Series["Series2"]->Points->Clear();

				 int a1, a2, a3, b1,b2,b3;
				 double k;
				 a = Convert::ToDouble(textBox12->Text);
				 l = Convert::ToDouble(textBox14->Text);
				 T = Convert::ToDouble(textBox13->Text);
				 h = Convert::ToDouble(textBox15->Text);
				 tau = Convert::ToDouble(textBox16->Text);

				 a1 = Convert::ToInt32(textBox1->Text);
				 a2 = Convert::ToInt32(textBox2->Text);
				 a3 = Convert::ToInt32(textBox3->Text);
				 b1 = Convert::ToInt32(textBox4->Text);
				 b2 = Convert::ToInt32(textBox5->Text);
				 b3 = Convert::ToInt32(textBox6->Text);

				 B1 = Convert::ToInt32(textBox7->Text);
				 B2 = Convert::ToInt32(textBox8->Text);
				 B3 = Convert::ToInt32(textBox9->Text);
				 B4 = Convert::ToInt32(textBox10->Text);
				 B5 = Convert::ToInt32(textBox11->Text);
				 int countI = l / h;
				 int countJ = T / tau;

				 //создание сетки с узлами
				 y = new double*[countI];

	
				 for (int i = 0; i < countI; i++)
					 y[i] = new double[countJ];


				 //0 слой
				 for (int i = 0; i < countI; i++){
					 y[i][0] = phi(a1, a2, a3, i*h);
				 }

				 for (int i = 0; i < countI; i++){
					 chart1->Series["Series1"]->Points->AddXY(i, y[i][0]);
					// chart1->Series["Series2"]->Points->AddXY(i, y[i][countJ - 2]);
				 }
				 //1 слой
				 for (int i = 1; i < countI-1; i++){
					 y[i][1] = phi(a1, a2, a3, i*h) + psi(b1, b2, b3, i*h)*tau + tau*a*(phi(a1, a2, a3, (i - 1)*h) - 2 * phi(a1, a2, a3, i*h) +
						 phi(a1, a2, a3, (i + 1)*h)) / (2 * h*h) + y[i][0]*bx(i*h) + y[i][0]*Simpson(&bx,l);
				 }

				 y[0][1] = (4 * y[1][1] - y[2][1]) / 3;
				 y[countI - 1][1] = (4 * y[countI - 2][1] - y[countI - 3][1]) / 3;
				 
				 double **A = new double*[countI];
				 for (int i = 0; i < countI; i++){
					 A[i] = new double[countI];
				 }

				 double *b = new double[countI];

				 for (int j = 2; j < countJ - 1; j++){
					 A = GenerateMatrix(a, tau, h, countI);
					 b = GenerateVector(y, &bx, tau, h, countI, j);


					 double* gaussSolution = new double[countI];
					 gaussSolution = gauss(A, b, countI);


					 for (int j1 = 0; j1 < countI; j1++) {
						 y[j1][j] = gaussSolution[j1];
					 }

					
				 }


				 for (int i = 0; i < countI; i++){
					// chart1->Series["Series1"]->Points->AddXY(i, y[i][0]);
					 chart1->Series["Series2"]->Points->AddXY(i, y[i][countJ-2]);
				 }
	}
private: System::Void button2_Click(System::Object^  sender, System::EventArgs^  e) {
			 	 SaveFileDialog^ saveFileDialog1 = gcnew SaveFileDialog;
				 saveFileDialog1->Filter = "Image Files(*.bmp)|*.bmp|Image Files(*.jpg)|*.jpg|Image Files(*.gif)|*.gif|Image Files(*.png)|*.png|All files (*.*)|*.*";
				 saveFileDialog1->FilterIndex = 1;
				 saveFileDialog1->RestoreDirectory = true;
				 if (saveFileDialog1->ShowDialog() == System::Windows::Forms::DialogResult::OK)
				 {
					chart1->SaveImage(saveFileDialog1->FileName, System::Drawing::Imaging::ImageFormat::Bmp);
				 }
		 }
private: System::Void MyForm_Load(System::Object^  sender, System::EventArgs^  e) {
}
private: System::Void chart1_Click(System::Object^  sender, System::EventArgs^  e) {
}
private: System::Void textBox1_TextChanged(System::Object^  sender, System::EventArgs^  e) {
}
private: System::Void groupBox1_Enter(System::Object^  sender, System::EventArgs^  e) {
}
private: System::Void label21_Click(System::Object^  sender, System::EventArgs^  e) {
}
private: System::Void label13_Click(System::Object^  sender, System::EventArgs^  e) {
}
};
}
