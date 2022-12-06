using System;
using System.Collections.Generic;
using Optical;
using Shapes;

namespace MathsTest
{
    class Program
    {
        static void Main(string[] args)
        {

            /*I want to create a lens. A standard lens must have two spherical
            faces. A planar lens has one spherical mirror and a flat face.*/

            SellmeierCoefficients sellmeier = new SellmeierCoefficients(1.03961212, 0.00600069867, 0.231792344, 0.0200179144, 0, 103);

            Element element = new Element(new Point<double>(0, 0), 10, 20, new Spherical_Refractive(20), new Spherical_Refractive(-50), sellmeier);


            List<Ray> rays = new()
            {
                new Ray(0.765, 0, new Point<double>(-10, 8))
            };

            foreach (Ray ray in rays)
            {
                

                Console.WriteLine("Lines:\n");

                //foreach (Line line in ray.Lines)
                //{
                //    Console.WriteLine($"y={line.B}x + {line.C}");
                //}

                

            }

            Console.WriteLine("\nPolygon Points:\n");

            foreach (Point<double> point in element.Polygon)
            {
                Console.WriteLine($"{point.X}, {point.Y}");
            }

            Console.WriteLine();
        }
    }
}
