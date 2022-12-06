using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace Shapes
{
    

    public class Point<T>
    {
    
        public Point(T x, T y)
        {
            X = x;
            Y = y;
        }
    
        public T X { get; set; }
    
        public T Y { get; set; }

        static public Point<double> AngleToVector(double Angle)
        {
            return new Point<double>(Math.Cos(-Angle), Math.Sin(-Angle));
        }
    }
    
    public struct PointVector
    {
        public PointVector(Point<double> point, Point<double> vector)
        {
            Point = point;
            Vector = vector;
        }
    
        public Point<double> Point { get; set; }
        public Point<double> Vector { get; set; }

        /// <summary>
        /// Converts the PointVector to a Line.
        /// </summary>
        /// <param name="line"></param>
        /// <returns></returns>
        public Line ToLine()
        {
            double angle, b, c;

            angle = -Math.Atan2(this.Vector.Y, this.Vector.X);

            b = Math.Tan(angle);
            c = this.Point.Y - (this.Point.X * b);

            return new Line(b, c);
        }
    }
    
    /// <summary>
    /// A circle as defined by (x - a)² + (y - b)² = r².
    /// </summary>
    public struct Circle
    {
        public double Radius { get; set; }
    
        public Point<double> Position { get; set; }
    
        public Circle(double radius, Point<double> position)
        {
            Radius = radius;
            Position = position;
        }
    }
    
    public struct Line
    {
        public double B { get; set; }
        public double C { get; set; }
    
        public Line(double b, double c)
        {
            B = b;
            C = c;
        }

        /// <summary>
        /// Gets the angle component of a Line
        /// </summary>
        /// <param name="line"></param>
        /// <returns></returns>
        public double GetAngle(Line line) { return Math.Asin(line.B); }
        /// <summary>
        /// Converts an angle in Radians to a Vector.
        /// </summary>
        /// <param name="Angle"></param>
        /// <returns></returns>
    }

    

}
