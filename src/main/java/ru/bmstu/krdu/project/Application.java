package ru.bmstu.krdu.project;

import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import ru.bmstu.krdu.project.dto.AstraResult;
import ru.bmstu.krdu.project.fileUtil.AstraReader;
import ru.bmstu.krdu.project.fileUtil.AstraUtil;

import java.text.DecimalFormat;
import java.util.*;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Function;

import ru.bmstu.krdu.project.dto.Vector;

public class Application {

    private static final double FOTH_FOR_OPTIMAL = 189;
    private static final double[] ALPHAS = new double[]{0.85, 0.86, 0.87, 0.88, 0.89, 0.9};
    private static final double GG_ALPHA = 13;
    private static final double[] ALPHAS_PR = new double[]{0.375, 0.38, 0.385};
    private static final double KM0 = 3.40496;
    private static final double L_PRIW = 1.7;

    private static final double D_KS = 0.20;
    private static final double VOLUME_UP_TO_CR = 2810000 * Math.pow(10, -9);
    private static final double P = 83360; //тяга, Ньютоны
    private static final double MYA_DIV_MSUM = 0.95; //отношение расхода компонентов ядра к общему
    private static final double MPR_DIV_MSUM = 0.05; //отношение расхода компонентов пристенка к общему

    private static final DecimalFormat decimal0Format = new DecimalFormat("0");
    private static final DecimalFormat decimal1Format = new DecimalFormat("0.0");
    private static final DecimalFormat decimal2Format = new DecimalFormat("0.00");
    private static final DecimalFormat decimal3Format = new DecimalFormat("0.000");
    private static final DecimalFormat decimal4Format = new DecimalFormat("0.0000");
    private static final DecimalFormat decimal5Format = new DecimalFormat("0.00000");

    private static final double BETA_A = 12 * Math.PI / 180; // рад
    private static final double BETA_M = 0.7178; // рад
    private static final double X_A = 33.92;

    public static void main(String[] args) {
        AstraReader astraReader = new AstraReader();

        var resInKS = astraReader.readParams("optimalAlpha\\KRDU.RES");
        var sectionResult = astraReader.readParams("section\\KRDU.RES");
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
        double r1 = 1.2 * dKs / 2;
        double r2 = 1.2 * dCr / 2;
        System.out.printf("R1 = 1.4 * Rk = 1.4 * %s / 2 = %s м\n", decimal3Format.format(dKs), decimal3Format.format(r1));
        System.out.printf("R2 = 1.3 * Rkr = 1.3 * %s / 2 = %s м\n", decimal3Format.format(dCr), decimal3Format.format(r2));

        var funcsOfNozzle = findNozzle(r1, r2, dKs / 2, dCr / 2, aCoefOfParab, bCoefOfParab, lRas, VOLUME_UP_TO_CR, true);
        printPoints(funcsOfNozzle, 200);
        var sectionPoints = findXForAstraRes(funcsOfNozzle, sectionResult, mFlowTotal);

        double phiParallel = (1 + Math.cos(BETA_A)) / 2;
        System.out.printf("phiP = (1 + cos(%s)) / 2 = %s\n", Math.toDegrees(BETA_A), decimal3Format.format(phiParallel));
        FrictionResult frictionResult = getFrictionForces(sectionPoints);

        double totalFrictionForce = Arrays.stream(frictionResult.getFrictionForces()).sum();
        System.out.printf("Ptr = sum(diffPtr) = %s Н\n", decimal1Format.format(totalFrictionForce));

        double phiFriction = (P - totalFrictionForce) / P;
        System.out.printf("phiTr = (P - Ptr) / P = (%s - %s) / %s = %s\n", P, decimal1Format.format(totalFrictionForce), P, decimal3Format.format(phiFriction));

        double phiNozzle = phiFriction * phiParallel;
        System.out.printf("phiS = phiP * phiTr = %s * %s = %s\n", decimal3Format.format(phiParallel), decimal3Format.format(phiFriction), decimal3Format.format(phiNozzle));

        double mFlowInCenter = MYA_DIV_MSUM * mFlowTotal;
        double mFlowOnWall = MPR_DIV_MSUM * mFlowTotal;
        System.out.printf("mtya = %s * msum = %s * %s = %s кг/c\n", MYA_DIV_MSUM, MYA_DIV_MSUM, decimal1Format.format(mFlowTotal), decimal1Format.format(mFlowInCenter));
        System.out.printf("mtpr = %s * msum = %s * %s = %s кг/c\n", MPR_DIV_MSUM, MPR_DIV_MSUM, decimal1Format.format(mFlowTotal), decimal1Format.format(mFlowOnWall));

        AstraResult ggAstraRes = astraReader.readParams("gg\\KRDU.RES").stream().findFirst().get();
        System.out.printf("alphaGg = %s\n", GG_ALPHA);
        System.out.printf("Tgg = %s K\n", decimal1Format.format(ggAstraRes.getT()));

        double kmInCenter = optMap.get("alpha") * KM0;
        double mFlowFuelCenter = mFlowInCenter / (1 + KM0);
        System.out.printf("mGYa = mtya / (1 + kmya) = %s / (1 + %s) = %s кг/c\n", decimal2Format.format(mFlowInCenter) , decimal2Format.format(kmInCenter), decimal2Format.format(mFlowFuelCenter));

        double mFlowOxInCenter = mFlowInCenter - mFlowFuelCenter;
        System.out.printf("mOkYa = mtya - mGYa = %s - %s = %s кг/c\n", decimal2Format.format(mFlowInCenter), decimal2Format.format(mFlowFuelCenter), decimal2Format.format(mFlowOxInCenter));

        double mFlowFuelGGCenter = mFlowOxInCenter / (GG_ALPHA * KM0);
        System.out.printf("mGGGYa = mOkYa / KmGG = %s / %s = %s кг/с\n", decimal2Format.format(mFlowOxInCenter), decimal2Format.format(GG_ALPHA * KM0), decimal2Format.format(mFlowFuelGGCenter));

        double mFlowFuelFluidCenter = mFlowFuelCenter - mFlowFuelGGCenter;
        System.out.printf("mGJidYa = mGYa - mGGGYa = %s - %s = %s кг/с\n", decimal2Format.format(mFlowFuelCenter), decimal2Format.format(mFlowFuelGGCenter), decimal2Format.format(mFlowFuelFluidCenter));

        double mFlowGGYa = mFlowFuelGGCenter + mFlowOxInCenter;
        System.out.printf("mGGYa = mGGGYa + mOkYa = %s + %s = %s кг/с\n", decimal2Format.format(mFlowFuelGGCenter), decimal2Format.format(mFlowOxInCenter), decimal2Format.format(mFlowGGYa));
    }

    private static FrictionResult getFrictionForces(List<Vector<AstraResult>> sectionPoints) {
        double[] squares = new double[sectionPoints.size() - 1];
        double[] cCoef = new double[sectionPoints.size() - 1];
        double[] betas = new double[sectionPoints.size() - 1];
        double[] diffX = new double[sectionPoints.size() - 1];
        double[] avgR = new double[sectionPoints.size() - 1];
        double[] m = new double[sectionPoints.size() - 1];
        double[] w = new double[sectionPoints.size() - 1];
        double[] rho = new double[sectionPoints.size() - 1];
        double[] frictionForces = new double[sectionPoints.size() - 1];
        double cFDefault = 0.006;
        double r = 0.89;
        System.out.println("------------------------------------------------------------------------");
        System.out.println("beta_ai - Rcr - diffXi - diffSi - Mi - Cfi - Wi - rhoi - diffPi");
        for (int i = 0; i < squares.length; i++) {
            Vector<AstraResult> prevVector = sectionPoints.get(i);
            Vector<AstraResult> currVector = sectionPoints.get(i + 1);
            diffX[i] = (currVector.getX() - prevVector.getX());
            avgR[i] = (currVector.getY().getD() + prevVector.getY().getD()) / 4;
            double beta = Math.atan(((currVector.getY().getD() - prevVector.getY().getD()) / 2) / diffX[i]);
            if(i == squares.length - 1) {
                beta = BETA_A;
            }
            betas[i] = Math.abs(Math.min(beta, BETA_M));

            squares[i] = 2 * Math.PI * diffX[i] * avgR[i];
            m[i] = currVector.getY().getM();
            w[i] = currVector.getY().getW();
            rho[i] = currVector.getY().getRho();
            cCoef[i] = cFDefault * Math.pow(1 + r * (((currVector.getY().getK() - 1)) / 2) * Math.pow(m[i], 2), -0.55);
            frictionForces[i] = squares[i] * Math.cos(betas[i]) * cCoef[i] * (rho[i] * Math.pow(w[i], 2)) / 2;
            System.out.printf("%-8s %-6s %-6s %-6s %-6s %-6s %-6s %-6s %-6s\n",
                    decimal2Format.format(Math.toDegrees(betas[i])),
                    decimal3Format.format(avgR[i]),
                    decimal4Format.format(diffX[i]),
                    decimal4Format.format(squares[i]),
                    decimal3Format.format(m[i]),
                    decimal4Format.format(cCoef[i]),
                    decimal1Format.format(w[i]),
                    decimal2Format.format(rho[i]),
                    decimal2Format.format(frictionForces[i])
            );
        }
        System.out.println("------------------------------------------------------------------------");
        return new FrictionResult(squares, cCoef, betas, diffX, avgR, m, w, rho, frictionForces);
    }

    private static List<Vector<AstraResult>> findXForAstraRes(Map<Double, Function<Double, Double>> funcsOfNozzle, List<AstraResult> sectionOriginal, double mSumFlow) {

        List<AstraResult> section = sectionOriginal.subList(1, sectionOriginal.size() - 1);

        SplineInterpolator interpolator = new SplineInterpolator();
        sectionOriginal.forEach(ar -> {
            ar.setF(ar.getFflow() * mSumFlow);
            ar.setD(Math.sqrt(4 * ar.getF() / Math.PI));
            ar.setRho((ar.getP() * Math.pow(10, 6)) / (ar.getR() * ar.getT()));
        });
        int idxOfCritic = 0;
        for (int i = 1; i < section.size(); i++) {
            if (section.get(i).getFoth() == 1) {
                idxOfCritic = i;
                break;
            }
        }
        List<AstraResult> beforeCritic = section.subList(0, idxOfCritic + 1);
        List<AstraResult> afterCritic = section.subList(idxOfCritic + 1, section.size());

        double[] boarder = new double[5];
        AtomicInteger j = new AtomicInteger(1);
        funcsOfNozzle.keySet().forEach(d -> {
            boarder[j.get()] = d;
            j.getAndIncrement();
        });

        double[] xBeforeCritic = new double[60];
        double[] yBeforeCritic = new double[60];
        double[] xAfterCritic = new double[60];
        double[] yAfterCritic = new double[60];

        for (int i = 1; i < boarder.length; i++) {
            if (i <= 3) {
                double step = (boarder[i] - boarder[i - 1]) / 20;
                double currX = boarder[i - 1];
                var func = funcsOfNozzle.get(boarder[i]);
                for (int k = (i - 1) * 20; k < i * 20; k++) {
                    xBeforeCritic[k] = currX;
                    yBeforeCritic[k] = (-1) * func.apply(currX) + 0.0001; //костыль
                    currX += step;
                }
            } else {
                double step = (boarder[i] - boarder[i - 1]) / 60;
                double currX = boarder[i - 1];
                var func = funcsOfNozzle.get(boarder[i]);
                for (int k = 0; k < yAfterCritic.length; k++) {
                    xAfterCritic[k] = currX;
                    yAfterCritic[k] = func.apply(currX) + 0.10 * k / yAfterCritic.length; //костыль
                    currX += step;
                }
                break;
            }
        }
        /*System.out.println(Arrays.stream(beforeCritic.stream().map(a -> a.getD() / 2).toArray()).toList());
        System.out.println(Arrays.stream(afterCritic.stream().map(a -> a.getD() / 2).toArray()).toList());
        System.out.println(Arrays.toString(yBeforeCritic));*/
        var beforeCrFuncRadiusToX = interpolator.interpolate(yBeforeCritic, xBeforeCritic);
        var afterCrFuncRadiusToX = interpolator.interpolate(yAfterCritic, xAfterCritic);

        List<Vector<AstraResult>> result = new LinkedList<>();
        beforeCritic.forEach(a -> result.add(new Vector<>(beforeCrFuncRadiusToX.value((-1) * a.getD() / 2), a)));
        afterCritic.forEach(a -> result.add(new Vector<>(afterCrFuncRadiusToX.value(a.getD() / 2), a)));
        result.add(0, new Vector<>(0, sectionOriginal.get(0)));
        result.add(result.size(), new Vector<>(boarder[boarder.length - 1], sectionOriginal.get(sectionOriginal.size() - 1)));
        return result;
    }

    private static Map<Double, Function<Double, Double>> findNozzle(double r1, double r2, double rKs, double rCr, double aCoef, double bCoef, double lRas, double volume, boolean printIterations) {
        double epsilon = 1000;
        SimpsonIntegrator integrator = new SimpsonIntegrator();
        if (printIterations) {
            System.out.println("----------------------------------");
        }
        double start = Math.sqrt((r1 + r2) * (r1 + r2) - (rKs - r1 - rCr - r2) * (rKs - r1 - rCr - r2));
        for (int i = (int) ((start) * 1000) + 1; i < 1000; i++) {
            var func = nozzle(r1, r2, (double) (i * 0.001), rKs, rCr, aCoef, bCoef, lRas);
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
            if (printIterations) {
                System.out.printf("i = %s мм, vol = %s м^3\n", i, decimal5Format.format(vol));
            }
            if (Math.abs(vol - volume) < 0.0001) {
                System.out.println("lSuj = " + (double) (i * 0.001) + " м");
                if (printIterations) {
                    System.out.println("----------------------------------");
                }
                return func;
            }
        }
        return null;
    }

    private static Map<Double, Function<Double, Double>> nozzle(double r1, double r2, double lSuj, double rKs, double rCr, double aCoef, double bCoef, double lRas) {
        Map<Double, Function<Double, Double>> nozzle = new LinkedHashMap<>();
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

        nozzle.put(x0, x -> Math.sqrt(r1 * r1 - x * x) + rKs - r1);
        nozzle.put(x1, x -> k * x + b);
        nozzle.put(lSuj, x -> (-1) * Math.sqrt(r2 * r2 - (x - lSuj) * (x - lSuj)) + rCr + r2);
        double cCoef = lSuj - aCoef * rCr * rCr - bCoef * rCr;

        nozzle.put(lSuj + lRas, x -> ((-1) * bCoef + Math.sqrt(bCoef * bCoef - 4 * aCoef * (cCoef - x))) / (2 * aCoef));
        return nozzle;
    }

    private static List<Vector<Double>> printPoints(Map<Double, Function<Double, Double>> functionMap, int numOfPoints) {
        System.out.println("----------------------");
        System.out.println("СОПЛО БЕЗ КС");
        List<Vector<Double>> vectorList = new ArrayList<>();
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
            vectorList.add(new Vector<>(x, y));
            System.out.printf("%s  %s\n", x * 1000, y * 1000);
        }
        System.out.println("----------------------");
        return vectorList;
    }

    private static class FrictionResult {
        private double[] squares;
        private double[] cCoef;
        private double[] betas;
        private double[] diffX;
        private double[] avgR;
        private double[] m;
        private double[] w;
        private double[] rho;
        private double[] frictionForces;

        public FrictionResult(double[] squares, double[] cCoef, double[] betas, double[] diffX, double[] avgR, double[] m, double[] w, double[] rho, double[] frictionForces) {
            this.squares = squares;
            this.cCoef = cCoef;
            this.betas = betas;
            this.diffX = diffX;
            this.avgR = avgR;
            this.m = m;
            this.w = w;
            this.rho = rho;
            this.frictionForces = frictionForces;
        }

        public double[] getSquares() {
            return squares;
        }

        public double[] getcCoef() {
            return cCoef;
        }

        public double[] getBetas() {
            return betas;
        }

        public double[] getDiffX() {
            return diffX;
        }

        public double[] getAvgR() {
            return avgR;
        }

        public double[] getM() {
            return m;
        }

        public double[] getW() {
            return w;
        }

        public double[] getRho() {
            return rho;
        }

        public double[] getFrictionForces() {
            return frictionForces;
        }
    }

}
