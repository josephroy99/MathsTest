using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Reflection;

namespace OpticalOldNope
{
    /// <summary>
    /// An Element contains two surfaces and material coefficients.
    /// </summary>
    class Element
    {
        public List<Point<double>> Polygon { get; set; }

        SellmeierCoefficients c;

        /// <summary>
        /// Element Position
        /// </summary>
        public Point<double> Position { get; set; }

        double _depth;
        /// <summary>
        /// Total Element depth at the centre.
        /// </summary>
        public double Depth { get; }

        /// <summary>
        /// Surface property
        /// </summary>
        public Surface[] Surfaces { get; }

        /// <summary>
        /// Lens diameter.
        /// </summary>
        public double Diameter { get; }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="Position">Element Position</param>
        /// <param name="Depth">Total Element depth at the centre.</param>
        /// <param name="Surface1"></param>
        /// <param name="Surface2"></param>
        /// <param name="c">Sellmeier Coefficients</param>
        public Element(Point<double> Element_Position, double Element_Depth, double Element_Diameter, Surface Surface1, Surface Surface2, SellmeierCoefficients C)
        {
            Polygon = new List<Point<double>>();
            Surfaces = new Surface[2];
            Surfaces[0] = Surface1;
            Surfaces[0].Parent = this;
            Surfaces[0].isFront = true;
            Surfaces[1] = Surface2;
            Surfaces[1].Parent = this;
            Surfaces[1].isFront = false;

            Depth = Element_Depth;
            Diameter = Element_Diameter;
            Position = Element_Position;
            c = C;
        }

        /// <summary>
        /// Takes a Ray as a parameter and calculates the surface hit points through the element.
        /// </summary>
        /// <param name="ray"></param>
        public void GeneratePolygon()
        {
            foreach (Surface surface in Surfaces)
            {
                surface.UpdateSurfacePositions();
                surface.GeneratePolyLine();
                Polygon.AddRange(surface.PolyLine);
            }
        }

        public void Flip()
        {
            var temp = Surfaces[0];
            Surfaces[0] = Surfaces[1];
            Surfaces[1] = temp;

            Surfaces[0].isFront = true;
            Surfaces[1].isFront = false;

        }
    }


    abstract class Surface
    {
        public Surface()
        {
            PolyLine = new List<Point<double>>();
            sign = isFront ? 1 : -1;
        }

        public bool isFront { get; set; }

        //This is the sign for front or rear. 1 if it is the front, -1 if it is the rear.
        protected double sign;

        /// <summary>
        /// Parent Element of the Surface.
        /// </summary>
        public Element Parent { get; set; }
        /// <summary>
        /// All points go clockwise from the element centre.
        /// </summary>
        public List<Point<double>> PolyLine { get; }
        public QuadraticFormula QuadraticFormula { get; }

        /// <summary>
        /// This function calculates the lines and points through a reflective element and applies this to the Ray given the glass properties.
        /// </summary>
        /// <param name="parent">Element that contains the surface</param>
        /// <param name="sign">This defines if it is the front or rear surface, -1 if it is the front, 1 if it is the rear.</param>
        /// <param name="ray">The Ray that will contain the results</param>
        public virtual void Calculate(int sign, Ray ray) { }
        /// <summary>
        /// This function calculates the lines and points through a refractive element and applies this to the Ray given the glass properties.
        /// </summary>
        /// <param name="parent">Element that contains the surface</param>
        /// <param name="sign">This defines if it is the front or rear surface, -1 if it is the front, 1 if it is the rear.</param>
        /// <param name="ray">This is the Ray that will contain the results</param>
        /// <param name="c">Sellmeier Coefficients</param>
        public virtual void Calculate(Point<double> hitPosition, Ray ray, SellmeierCoefficients s) { }

        //Functions
        /// <summary>
        /// Converts a PointVector to a Line.
        /// </summary>
        /// <param name="line"></param>
        /// <returns></returns>
        public Line PointVectorToLine(PointVector pointVector)
        {
            double angle, b, c;

            angle = -Math.Atan2(pointVector.Vector.Y, pointVector.Vector.X);

            b = Math.Tan(angle);
            c = pointVector.Point.Y - (pointVector.Point.X * b);

            return new Line(b, c);
        }
        /// <summary>
        /// Gets the angle component of a Line
        /// </summary>
        /// <param name="line"></param>
        /// <returns></returns>
        protected double GetAngleOfLine(Line line) { return Math.Asin(line.B); }
        /// <summary>
        /// Converts an angle in Radians to a Vector.
        /// </summary>
        /// <param name="Angle"></param>
        /// <returns></returns>
        protected Point<double> AngleToVector(double Angle)
        {
            return new Point<double>(Math.Cos(-Angle), Math.Sin(-Angle));
        }
        public virtual QuadraticFormula GetQuadraticFormula(PointVector pointVector)
        {
            return new QuadraticFormula();
        }

        public abstract Point<double> FindIntersection(PointVector pointVector);

        public abstract void UpdateSurfacePositions();


        public abstract void GeneratePolyLine();
    }

    /// <summary>
    /// Contains the function to find intersection of two Lines.
    /// </summary>
    abstract class Flat_Surface : Surface
    {
        public override void UpdateSurfacePositions()
        {
            sign = isFront ? 1 : -1;
            surface = new Line(0, -sign * Parent.Depth / 2);
        }

        //This line is defined as x = by+c. It will need converting to y = bx+c.
        protected Line surface;

        /// <summary>
        /// This function takes a PointVector and checks if there a collision between itself and a Line. <Unchecked>
        /// </summary>
        /// <param name="pointVector">A PointVector where (1,0) is 0 rad.</param>
        /// <returns></returns>
        public override Point<double> FindIntersection(PointVector pointVector)
        {
            //Converts point vector to a line.
            Line line = PointVectorToLine(pointVector);
            //This is where the conversion takes place to make it y = bx + c.
            Line surfaceLine = new Line(1 / surface.B, (1 / surface.B) * surface.C);

            double x, y, t, u;

            //If C is infinity it means the line is vertical and can be solved by inserting x into the equation.
            if (Double.IsInfinity(surfaceLine.C))
            {
                y = (line.B * surface.C) + line.C;
                x = surface.C;
            }
            else //If not it means it is angled and therefore needs solving by substitution.
            {
                t = line.B + surfaceLine.B;
                u = line.C + surfaceLine.C;

                x = t / u;
                y = (line.B * x) + line.C;
            }

            return new Point<double>(x, y);
        }

        /// <summary>
        /// Generates a PolyLine, which can be merged with another to make a full element shape.
        /// </summary>
        public override void GeneratePolyLine()
        {
            PolyLine.Add(new Point<double>(-sign * Parent.Depth / 2, sign * -Parent.Diameter / 2));
            PolyLine.Add(new Point<double>(-sign * Parent.Depth / 2, sign * Parent.Diameter / 2));
        }
    }

    abstract class Spherical_Surface : Surface
    {
        public override void UpdateSurfacePositions()
        {
            sign = (isFront ? 1 : -1);
            surface.Position = new Point<double>(sign * (surface.Radius - (Parent.Depth / 2)), 0);
        }

        protected Circle surface;

        public override Point<double> FindIntersection(PointVector pointVector)
        {
            Line line = PointVectorToLine(pointVector);
            QuadraticFormula f = GetQuadraticFormula(pointVector);
            double x, y;
            x = f.X;
            y = (x * line.B) + line.C;

            return new Point<double>(x, y);
        }

        public override QuadraticFormula GetQuadraticFormula(PointVector pointVector)
        {
            //Convert point vector to a line.
            Line line = PointVectorToLine(pointVector);

            double xs2 = (surface.Position.X * -2);
            double ys2 = (surface.Position.Y * -2);

            double a, b, c;

            int sign = surface.Position.X > 0 ? -1 : 1;

            //Creates a quadratic formuala via subsition of a Line and a Circle equation.
            a = 1 + (Math.Pow(line.B, 2));
            b = xs2 + (ys2 * line.B) + (line.B * line.C * 2);
            c = Math.Pow(line.C, 2) + (line.C * ys2) + Math.Pow(surface.Position.X * -1, 2) + Math.Pow(surface.Position.Y * -1, 2) - Math.Pow(surface.Radius, 2);

            //Solves the quadratic formula via a quadratic equation.

            return new QuadraticFormula(a, b, c);
        }

        /// <summary>
        /// Generates a PolyLine, which can be merged with another to make a full element shape.
        /// </summary>
        public override void GeneratePolyLine()
        {
            Point<double> arcCentre;
            double arcTotalSweep, arcResolutionSweep;
            int resolution = 30;

            arcTotalSweep = Math.Asin((Parent.Diameter / 2) / surface.Radius) * 2;
            arcResolutionSweep = arcTotalSweep / (resolution - 1);

            arcCentre = new Point<double>(surface.Position.X, surface.Position.Y);

            for (int i = 0; i < resolution; i++)
            {
                double currentAngle = (sign * arcTotalSweep / 2) + -sign * (arcResolutionSweep * i);
                double x, y;

                x = Math.Cos(currentAngle) * -sign * surface.Radius + arcCentre.X;
                y = Math.Sin(currentAngle) * surface.Radius;


                PolyLine.Add(new Point<double>(x, y));
            }
        }
    }

    class Flat_Refractive : Flat_Surface
    {
        public Flat_Refractive()
        {

        }

        public override void Calculate(Point<double> hitPosition, Ray ray, SellmeierCoefficients s)
        {

        }
    }

    class Spherical_Refractive : Spherical_Surface
    {
        public Spherical_Refractive(double Radius)
        {
            surface.Radius = Radius;
        }

        public override void Calculate(Point<double> hitPosition, Ray ray, SellmeierCoefficients s)
        {
            double n1, n2;
            if (isFront)
            {
                n1 = s.AirRefractiveIndex(ray.Wavelength); ;
                n2 = s.CalculateRefractiveIndex(ray.Wavelength);
            }
            else if (!isFront)
            {
                n1 = s.CalculateRefractiveIndex(ray.Wavelength);
                n2 = s.AirRefractiveIndex(ray.Wavelength); ;
            }
            else
            {
                n1 = 0;
                n2 = 0;
            }

            double lineAngle, angleOfIncidence, angleOfRefraction, angleOfRefractionFromCentreline, normalAngle;

            //need to convert vector to angle correctly
            int count = ray.PointVectors.Count - 1;
            lineAngle = -Math.Atan2(ray.PointVectors[count].Vector.Y, ray.PointVectors[count].Vector.X);
            normalAngle = Math.Asin(hitPosition.Y / (sign * surface.Radius));
            angleOfIncidence = normalAngle + lineAngle;
            angleOfRefraction = Math.Asin((n1 * Math.Sin(angleOfIncidence)) / n2);
            angleOfRefractionFromCentreline = angleOfRefraction - normalAngle;

            double b, c;
            b = Math.Tan(angleOfRefractionFromCentreline);
            c = hitPosition.Y - (hitPosition.X * b);

            ray.PointVectors.Add(new PointVector(new Point<double>(hitPosition.X, hitPosition.Y), AngleToVector(angleOfRefractionFromCentreline)));
            ray.Lines.Add(new Line(b, c));
        }
    }

    struct QuadraticFormula
    {
        public QuadraticFormula(double A, double B, double C)
        {
            a = A;
            b = B;
            c = C;
        }

        public double a { get; }
        public double b { get; }
        public double c { get; }

        public double X
        {
            get
            {
                double xPlus, xMinus;
                xPlus  = (-b + Math.Sqrt(Math.Pow(b, 2) - 4 * a * c)) / (2 * a);
                xMinus = (-b - Math.Sqrt(Math.Pow(b, 2) - 4 * a * c)) / (2 * a);

                return xPlus;
            }
        }
    }

    struct Line
    {
        public double B { get; set; }
        public double C { get; set; }

        public Line(double b, double c)
        {
            B = b;
            C = c;
        }
    }

    struct PointVector
    {
        public PointVector(Point<double> point, Point<double> vector)
        {
            Point = point;
            Vector = vector;
        }

        public Point<double> Point { get; set; }

        public Point<double> Vector { get; set; }
    }

    /// <summary>
    /// A circle as defined by (x - a)² + (y - b)² = r².
    /// </summary>
    struct Circle
    {
        public double Radius { get; set; }

        public Point<double> Position { get; set; }

        public Circle(double radius, Point<double> position)
        {
            Radius = radius;
            Position = position;
        }
    }

    class Point<T>
    {

        public Point(T x, T y)
        {
            X = x;
            Y = y;
        }

        public T X { get; set; }

        public T Y { get; set; }
    }

    /// <summary>
    /// A Ray contains a Wavelength, a list of hit points and a list of lines.
    /// </summary>
    struct Ray
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="Wavelenth">Wavelength in µm</param>
        /// <param name="Angle">Angle in Radians</param>
        /// <param name="Start_Position"></param>
        public Ray(double wavelenth, double angle, Point<double> Start_Position)
        {
            Wavelength = wavelenth;
            Lines = new List<Line>();
            PointVectors = new List<PointVector>();

            Lines.Add(new Line(Math.Tan(-angle), Start_Position.Y - (Start_Position.X * Math.Tan(-angle))));
            PointVectors.Add(new PointVector(Start_Position, new Point<double>(Math.Cos(-angle), Math.Sin(-angle))));
        }

        public double Wavelength { get; }

        public List<Line> Lines { get; }

        public List<PointVector> PointVectors { get; set; }

        


    }

    struct ElementHit
    {
        public Element element;
        public Surface surface;
    }

    struct SellmeierCoefficients
    {
        double _B1, _B2, _B3, _C1, _C2, _C3;

        public SellmeierCoefficients(double B1, double C1, double B2, double C2, double B3, double C3)
        {
            _B1 = B1;
            _B2 = B2;
            _B3 = B3;
            _C1 = C1;
            _C2 = C2;
            _C3 = C3;
        }

        public double CalculateRefractiveIndex(double Wavelength)
        {
            double wavelenths2 = Math.Pow(Wavelength, 2);
            double A1 = _C1 != 0 ? _B1 / (1 - _C1 / wavelenths2) : 0;
            double A2 = _C2 != 0 ? _B2 / (1 - _C2 / wavelenths2) : 0;
            double A3 = _C3 != 0 ? _B3 / (1 - _C3 / wavelenths2) : 0;

            return Math.Sqrt(1 + A1 + A2 + A3);
        }

        public double AirRefractiveIndex(double Wavelength)
        {
            double wavelenths2 = Math.Pow(Wavelength, 2);
            double A1 = 238.0185 != 0 ? 0.05792105 / (1 - 238.0185 / wavelenths2) : 0;
            double A2 = 57.362 != 0 ? 0.00167917 / (1 - 57.362 / wavelenths2) : 0;

            return Math.Sqrt(1 + A1 + A2);
        }
    }
}
