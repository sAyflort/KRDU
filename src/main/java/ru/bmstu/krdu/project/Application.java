package ru.bmstu.krdu.project;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.integration.TrapezoidIntegrator;
import ru.bmstu.krdu.project.fileUtil.AstraReader;
import ru.bmstu.krdu.project.fileUtil.AstraUtil;

import java.awt.*;
import java.awt.geom.Point2D;
import java.text.DecimalFormat;
import java.util.*;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Function;

public class Application {
    private static final double FOTH_FOR_OPTIMAL = 189;
    private static final double[] ALPHAS = new double[]{0.85, 0.86, 0.87, 0.88, 0.89, 0.9};
    private static final double[] ALPHAS_PR = new double[]{0.375, 0.38, 0.385};
    private static final double KM0 = 3.40496;
    private static final double L_PRIW = 1.7;

    private static final double D_KS = 0.15;
    private static final double VOLUME_UP_TO_CR = 1510000 * Math.pow(10, -9);
    private static final double P = 83360; //тяга, Ньютоны
    private static final double MYA_DIV_MSUM = 0.95; //отношение расхода компонентов ядра к общему
    private static final double MPR_DIV_MSUM = 0.05; //отношение расхода компонентов пристенка к общему

    private static final DecimalFormat decimal0Format = new DecimalFormat("0.");
    private static final DecimalFormat decimal1Format = new DecimalFormat("0.0");
    private static final DecimalFormat decimal2Format = new DecimalFormat("0.00");
    private static final DecimalFormat decimal3Format = new DecimalFormat("0.000");
    private static final DecimalFormat decimal4Format = new DecimalFormat("0.0000");

    private static final double BETA_A = 12 * Math.PI / 180; // рад
    private static final double BETA_M = 0.7178; // рад
    private static final double X_A = 33.92;

    public static void main(String[] args) {
        AstraReader astraReader = new AstraReader();

        var resInKS = astraReader.readParams("optimalAlpha\\KRDU.RES");
        var manyResult = astraReader.readParams("secheniya\\KRDU.RES");
        var optMap = AstraUtil.tableOfOptimum(resInKS, FOTH_FOR_OPTIMAL, ALPHAS, KM0);
        double iOpt = optMap.get("i") * 9.81;
        Map<String, Double> mapOfParamFromPr = AstraUtil.tableOfPristenok(astraReader.readParams("optimalAlpha\\KRDU_PR.RES"), FOTH_FOR_OPTIMAL, ALPHAS_PR, KM0);
        double iPr = mapOfParamFromPr.get("i") * 9.81;

        double iTotal = MYA_DIV_MSUM * iOpt + MPR_DIV_MSUM * iPr;
        System.out.printf("Km0 = %s\n", decimal3Format.format(KM0));
        System.out.printf("Iysum = %s * %s + %s * %s = %s м/c\n", MYA_DIV_MSUM, decimal2Format.format(iOpt), MPR_DIV_MSUM, decimal2Format.format(iPr), decimal2Format.format(iTotal));

        double mFlowTotal = P / (iTotal);
        System.out.printf("msum =  %s / %s = %s кг/c\n", P, decimal2Format.format(iTotal), decimal3Format.format(mFlowTotal));

        double fFlowCr = optMap.get("fFlowCr");
        double fCr = mFlowTotal * fFlowCr;
        System.out.printf("Fcr = %s * %s * 10^-4 = %s * 10^-3 м^2\n", decimal3Format.format(mFlowTotal), fFlowCr * 10000, decimal2Format.format((fCr * 1000)));

        double dCr = Math.sqrt(4 * fCr / Math.PI);
        System.out.printf("dCr = (4 * %s / 3.14)^0.5 = %s м\n", decimal2Format.format((fCr * 1000)), decimal3Format.format(dCr));

        double fFlowOut = optMap.get("fFlowOut");
        double fOut = mFlowTotal * fFlowOut;
        System.out.printf("Fa = %s * %s * 10^-4 = %s * 10^-3 м^2\n", decimal3Format.format(mFlowTotal), decimal2Format.format(fFlowOut * 10000), decimal2Format.format((fOut * 1000)));

        double dOut = Math.sqrt(4 * fOut / Math.PI);
        System.out.printf("da = (4 * %s / 3.14)^0.5 = %s м\n", decimal2Format.format((fOut * 1000)), decimal3Format.format(dOut));

        double fKs = Math.PI * D_KS * D_KS / 4;
        System.out.printf("Fks = %s * 10^-3 м^2\n", decimal2Format.format(fKs * 1000));

        double dKs = Math.sqrt(4 * fKs / Math.PI);
        System.out.printf("dks = (4 * %s / 3.14)^0.5 = %s м\n", decimal2Format.format((fKs * 1000)), decimal3Format.format(dKs));

        double volumeUpToCrWithKS = L_PRIW * fCr;
        System.out.printf("Vks + Vsuj = Lpriw * Fcr = %s * %s * 10^-3 = %s м^3\n", L_PRIW, decimal2Format.format((fCr * 1000)), decimal4Format.format(volumeUpToCrWithKS));

        System.out.printf("Vsuj = %s м^3\n", decimal4Format.format(VOLUME_UP_TO_CR));

        double lKs = (volumeUpToCrWithKS - VOLUME_UP_TO_CR) / fKs;

        System.out.printf("Lks = Vks / Fks = (%s - %s) / (%s * 10^-3) = %s м\n", decimal4Format.format(volumeUpToCrWithKS), VOLUME_UP_TO_CR, decimal3Format.format(fKs * 1000), decimal3Format.format(lKs));

        double k = optMap.get("k");
        double pOut = optMap.get("pa");
        double pKs = optMap.get("p");

        /*double step = (pKs - pOut) / 99;
        double[] pForAstra = new double[100];
        pForAstra[0] = pKs;
        for (int i = 1; i < pForAstra.length; i++) {
            pForAstra[i] = pKs - i * step;
        }
        pForAstra[pForAstra.length - 1] = pOut;
        System.out.println(Arrays.toString(pForAstra));*/

        double dOtn = dOut / dCr;
        double pOtn = pKs / pOut;

        System.out.println("k = " + k);
        System.out.println("dOtn = " + dOtn);
        System.out.println("pks/pa = " + pOtn);
        System.out.println("betaA = " + Math.toDegrees(BETA_A));
        System.out.println("betaM = " + decimal0Format.format(Math.toDegrees(BETA_M)));
        System.out.println("Xa = " + X_A);

        double aCoefOfParab = (Math.tan(Math.PI / 2 - BETA_A) - Math.tan(Math.PI / 2 - BETA_M)) / (dOut - dCr);
        double bCoefOfParab = Math.tan(Math.PI / 2 - BETA_M) - aCoefOfParab * dCr;
        Function<Double, Double> parab = x -> aCoefOfParab * x * x;

        double lRas = /*X_A * dCr / 2*/ parab.apply(dOut / 2) - parab.apply(dCr / 2);
        System.out.printf("Lras = Xa * Rcr = %s * %s / 2 = %s м\n", X_A, decimal4Format.format(dCr), decimal3Format.format(lRas));
        double r1 = 1.4 * dKs / 2;
        double r2 = 1.3 * dCr / 2;
        System.out.printf("R1 = 1.4 * Rk = 1.4 * %s / 2 = %s м\n", decimal3Format.format(dKs), decimal3Format.format(r1));
        System.out.printf("R2 = 1.3 * Rkr = 1.3 * %s / 2 = %s м\n", decimal3Format.format(dCr), decimal3Format.format(r2));
        var soplo = findSoplo(r1, r2, dKs / 2, dCr / 2, aCoefOfParab, bCoefOfParab, lRas, VOLUME_UP_TO_CR);
        System.out.println("");
        var points = printPoints(soplo, 200);
    }

    private static Map<Double, Function<Double, Double>> findSoplo(double r1, double r2, double rKs, double rCr, double aCoef, double bCoef, double lRas, double volume) {
        double epsilon = 1000;
        SimpsonIntegrator integrator = new SimpsonIntegrator();

        double start = Math.sqrt((r1 + r2) * (r1 + r2) - (rKs - r1 - rCr - r2) * (rKs - r1 - rCr - r2));

        for (int i = (int) (start) + 100; i < 1000; i++) {
            var func = soplo(r1, r2, (double) (i * 0.001), rKs, rCr, aCoef, bCoef, lRas);
            double vol = 0;
            double[] boarder = new double[5];
            AtomicInteger j = new AtomicInteger(1);
            func.keySet().forEach(d -> {
                boarder[j.get()] = d;
                j.getAndIncrement();
            });
            for (int k = 1; k < boarder.length - 1; k++) {
                int finalK = k;
                vol += integrator.integrate(100, x -> Math.pow(func.get(boarder[finalK]).apply(x), 2), boarder[k - 1], boarder[k]) * Math.PI;
            }
            //System.out.println("i, vol=" + i + ", " + vol);
            if (Math.abs(vol - volume) < 0.0001) {
                System.out.println("lSuj = " + (double) (i * 0.001) + " м");
                return func;
            }
        }


        return null;
    }

    private static Map<Double, Function<Double, Double>> soplo(double r1, double r2, double lSuj, double rKs, double rCr, double aCoef, double bCoef, double lRas) {
        Map<Double, Function<Double, Double>> profileSopla = new LinkedHashMap<>();
        double cx = (-1) * (rKs - r1 - (rCr + r2));
        double c = Math.sqrt(cx * cx + lSuj * lSuj);
        double l1 = (r1 * c / r2) / (1 + r1 / r2);
        double acosR1DivL1 = Math.acos(r1 / l1);
        double atanCxDivC = Math.atan(cx / lSuj);
        double arg = acosR1DivL1 + atanCxDivC;
        double x0 = r1 * Math.cos(arg);

        double k = (-1) * x0 / Math.sqrt(r1 * r1 - x0 * x0);
        double b = Math.sqrt(r1 * r1 - x0 * x0) - k * x0 + rKs - r1;

        double x1 = lSuj - Math.abs(k * r2) / Math.sqrt(1 + k * k);

        profileSopla.put(x0, x -> Math.sqrt(r1 * r1 - x * x) + rKs - r1);
        profileSopla.put(x1, x -> k * x + b);
        profileSopla.put(lSuj, x -> (-1) * Math.sqrt(r2 * r2 - (x - lSuj) * (x - lSuj)) + rCr + r2);
        double cCoef = lSuj - aCoef * rCr * rCr - bCoef * rCr;

        profileSopla.put(lSuj + lRas, x -> ((-1) * bCoef + Math.sqrt(bCoef * bCoef - 4 * aCoef * (cCoef - x))) / (2 * aCoef));
        return profileSopla;
    }

    private static List<Vector> printPoints(Map<Double, Function<Double, Double>> functionMap, int numOfPoints) {
        System.out.println("----------------------");
        System.out.println("СОПЛО БЕЗ КС");
        List<Vector> vectorList = new ArrayList<>();
        double[] boarder = new double[4];
        AtomicInteger j = new AtomicInteger(0);
        functionMap.keySet().forEach(d -> {
            boarder[j.get()] = d;
            j.getAndIncrement();
        });
        double step = boarder[boarder.length - 1] / numOfPoints;
        double x;
        int idxBoarder = 0;
        Function<Double, Double> func = functionMap.get(boarder[idxBoarder]);
        for (int i = 0; i <= numOfPoints; i++) {
            x = i * step;
            if (x > boarder[idxBoarder]) {
                //System.out.println("--------------------------------------------------");
                idxBoarder++;
                func = functionMap.get(boarder[idxBoarder]);
            }
            double y = func.apply(x);
            vectorList.add(new Vector(x, y));
            System.out.printf("%s  %s\n", x * 1000, y * 1000);
        }
        System.out.println("----------------------");
        return vectorList;
    }

    private static class Vector {
        private double x;
        private double y;

        public Vector(double x, double y) {
            this.x = x;
            this.y = y;
        }

        public double getX() {
            return x;
        }

        public void setX(double x) {
            this.x = x;
        }

        public double getY() {
            return y;
        }

        public void setY(double y) {
            this.y = y;
        }
    }
}
