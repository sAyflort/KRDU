package ru.bmstu.krdu.project;

import org.apache.commons.math3.analysis.integration.SimpsonIntegrator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.solvers.BrentSolver;
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

    private static final double RHO_OX = 1140;
    private static final double RHO_FUEL = 835;
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

    private static final Function<Double, Double> FRENEL_UNDER_INTEGRATE_FUNC = z -> Math.pow(Math.E, (-1) * z * z);

    private static final SimpsonIntegrator INTEGRATOR = new SimpsonIntegrator();
    private static final double BETA_A = 12 * Math.PI / 180; // рад
    private static final double BETA_M = 0.7178; // рад
    private static final double X_A = 33.92;

    private static int NUM_FUEL_SLICE = 6;
    private static int NUM_PR_SLICE = 2;
    private static int NUM_NOZZLE_PER_PR_SLICE = 48; // : 4
    private static double COEF_STEP_PR = 6.15; // > 1

    private static double DELTA_P_NOZZLE_PR = 1; // < DELTA_P_MAG_FUEL1
    private static double DELTA_P_NOZZLE_FUEL = 1; // < DELTA_P_MAG_FUEL1
    private static double DELTA_P_NOZZLE_GG = 1; // < DELTA_P_MAG_FUEL1
    private static double NU_NOZZLE_PR = 0.8; //0.75-0.85
    private static double NU_NOZZLE_GG = 0.8; //0.75-0.85
    private static double NU_NOZZLE_FUEL = 0.25; //0.1-0.5
    private static double TWO_ALPHA_NOZZLE_FUEL = Math.PI / 2;
    private static double PHI_NOZZLE_FUEL = 0.45;
    private static double A_NOZZLE_FUEL = 2.5;
    private static double R_INPUT_DIV_R_NOZZLE_NOZZLE_FUEL = 2.5;
    private static int NUM_OF_TAN_INPUT_NOZZLE_FUEL = 4;


    private static double P1_O = 1.2;
    private static double P1_G = 1.2;
    private static double ETA_OX = 0.6;
    private static double ETA_FUEL = 0.6;
    private static double ETA_T = 0.78;
    private static double DELTA_P_MAG_OX = 1.5;
    private static double DELTA_P_MAG_FUEL = 1.5;
    private static double DELTA_P_MAG_FUEL1 = 7;
    private static double DELTA_P_T1 = 1.5;
    private static double DELTA_P_T2 = 1;

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

        double lRas = parab.apply(dOut / 2) - parab.apply(dCr / 2);
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
        ggAstraRes.setRho((ggAstraRes.getP() * Math.pow(10, 6)) / (ggAstraRes.getR() * ggAstraRes.getT()));
        System.out.printf("alphaGg = %s\n", GG_ALPHA);
        System.out.printf("Tgg = %s K\n", decimal1Format.format(ggAstraRes.getT()));

        double kmInCenter = optMap.get("alpha") * KM0;
        double mFlowFuelCenter = mFlowInCenter / (1 + kmInCenter);
        System.out.printf("mGYa = mtya / (1 + kmya) = %s / (1 + %s) = %s кг/c\n", decimal2Format.format(mFlowInCenter), decimal2Format.format(kmInCenter), decimal2Format.format(mFlowFuelCenter));

        double mFlowOxInCenter = mFlowInCenter - mFlowFuelCenter;
        System.out.printf("mOkYa = mtya - mGYa = %s - %s = %s кг/c\n", decimal2Format.format(mFlowInCenter), decimal2Format.format(mFlowFuelCenter), decimal2Format.format(mFlowOxInCenter));

        double mFlowFuelGGCenter = mFlowOxInCenter / (GG_ALPHA * KM0);
        System.out.printf("mGGGYa = mOkYa / KmGG = %s / %s = %s кг/с\n", decimal2Format.format(mFlowOxInCenter), decimal2Format.format(GG_ALPHA * KM0), decimal2Format.format(mFlowFuelGGCenter));

        double mFlowFuelFluidCenter = mFlowFuelCenter - mFlowFuelGGCenter;
        System.out.printf("mGJidYa = mGYa - mGGGYa = %s - %s = %s кг/с\n", decimal2Format.format(mFlowFuelCenter), decimal2Format.format(mFlowFuelGGCenter), decimal2Format.format(mFlowFuelFluidCenter));

        double mFlowGGYa = mFlowFuelGGCenter + mFlowOxInCenter;
        System.out.printf("mGGYa = mGGGYa + mOkYa = %s + %s = %s кг/с\n", decimal2Format.format(mFlowFuelGGCenter), decimal2Format.format(mFlowOxInCenter), decimal2Format.format(mFlowGGYa));

        double pGG = getTNABalanceResult(ggAstraRes, mFlowOxInCenter, mFlowFuelCenter, pKs);

        double mFlowNozzlePr = mFlowOnWall / (NUM_NOZZLE_PER_PR_SLICE * NUM_PR_SLICE);
        System.out.printf("mFlowNozzlePr = mFlowOnWall / numOfNozzlePr = %s / %s = %s кг/c\n", decimal2Format.format(mFlowOnWall), NUM_NOZZLE_PER_PR_SLICE * NUM_PR_SLICE, decimal4Format.format(mFlowNozzlePr));


        double squareNozzlePr = mFlowNozzlePr / (NU_NOZZLE_PR * Math.sqrt(2 * DELTA_P_NOZZLE_PR * Math.pow(10, 6) * RHO_FUEL));
        System.out.printf("squareNozzlePr = mFlowNozzlePr / (NU_NOZZLE_PR * sqrt(2 * DELTA_P_NOZZLE_PR * Math.pow(10, 6) * RHO_FUEL)) = %s / (%s * sqrt(2 * %s * 1000000 * %s)) = %s * 10^-5 м^2\n",
                decimal5Format.format(mFlowNozzlePr), decimal2Format.format(NU_NOZZLE_PR), decimal2Format.format(DELTA_P_NOZZLE_PR), decimal2Format.format(RHO_FUEL), decimal3Format.format(squareNozzlePr * Math.pow(10, 5)));
        double dNozzlePr = Math.sqrt(4 * squareNozzlePr / Math.PI);
        System.out.printf("dNozzlePr = sqrt(4 * squareNozzlePr / PI) = sqrt(4 * %s * 10^-5 / PI) = %s мм\n", decimal3Format.format(squareNozzlePr * Math.pow(10, 5)), decimal2Format.format(dNozzlePr * 1000));

        //TODO нужна оптимизация
        MixedHead mixedHead = getMixedHead(NUM_PR_SLICE, NUM_NOZZLE_PER_PR_SLICE, NUM_FUEL_SLICE, dNozzlePr, dKs, true);
        double numOfNozzleFuel = mixedHead.getFuelNozzles().size();
        double numOfNozzleOx = mixedHead.getOxNozzles().size();
        System.out.printf("numOfNozzleOx = %s\n", numOfNozzleOx);
        System.out.printf("numOfNozzleFuel = %s\n", numOfNozzleFuel);

        double mFlowNozzleFuel = mFlowFuelFluidCenter / numOfNozzleFuel;
        System.out.printf("mFlowNozzleFuel = mFlowFuelFluidCenter / numOfNozzleFuel = %s / %s = %s кг/c\n", decimal2Format.format(mFlowFuelFluidCenter), decimal2Format.format(numOfNozzleFuel), decimal4Format.format(mFlowNozzleFuel));
        double mFlowNozzleGG = mFlowGGYa / numOfNozzleOx;
        System.out.printf("mFlowNozzleGG = mFlowGGYa / numOfNozzleOx = %s / %s = %s кг/c\n", decimal2Format.format(mFlowGGYa), decimal2Format.format(numOfNozzleOx), decimal4Format.format(mFlowNozzleGG));

        double squareNozzleFuel = mFlowNozzleFuel / (NU_NOZZLE_FUEL * Math.sqrt(2 * DELTA_P_NOZZLE_FUEL * Math.pow(10, 6) * RHO_FUEL));
        System.out.printf("squareNozzleFuel = mFlowNozzleFuel / (numOfNozzleFuel * sqrt(2 * DELTA_P_NOZZLE_FUEL * Math.pow(10, 6) * RHO_FUEL)) = %s / (%s * sqrt(2 * %s * 1000000 * %s)) = %s * 10^-5 м^2\n",
                decimal5Format.format(mFlowNozzleFuel), decimal2Format.format(numOfNozzleFuel), decimal2Format.format(DELTA_P_NOZZLE_FUEL), decimal2Format.format(RHO_FUEL), decimal3Format.format(squareNozzleFuel * Math.pow(10, 5)));
        double dNozzleFuel = Math.sqrt(4 * squareNozzleFuel / Math.PI);
        double rNozzleFuel = dNozzleFuel / 2;
        System.out.printf("dNozzleFuel = sqrt(4 * squareNozzleFuel / PI) = sqrt(4 * %s * 10^-5 / PI) = %s мм\n", decimal3Format.format(squareNozzleFuel * Math.pow(10, 5)), decimal2Format.format(dNozzleFuel * 1000));
        System.out.printf("Rvx / rc = %s\n", R_INPUT_DIV_R_NOZZLE_NOZZLE_FUEL);
        System.out.printf("A = %s\n", A_NOZZLE_FUEL);
        System.out.printf("phi = %s\n", PHI_NOZZLE_FUEL);

        double rInputNozzleFuel = R_INPUT_DIV_R_NOZZLE_NOZZLE_FUEL * (dNozzleFuel / 2);
        System.out.printf("rInputNozzleFuel = R_INPUT_DIV_R_NOZZLE_NOZZLE_FUEL * (dNozzleFuel / 2) = %s * (%s / 2) = %s мм\n", R_INPUT_DIV_R_NOZZLE_NOZZLE_FUEL, decimal2Format.format(dNozzleFuel * 1000), decimal2Format.format(rInputNozzleFuel * 1000));

        double rTanInputNozzleFuel = Math.sqrt((rInputNozzleFuel * rNozzleFuel) / (NUM_OF_TAN_INPUT_NOZZLE_FUEL * A_NOZZLE_FUEL));
        System.out.printf("rTanInput = sqrt((rInputNozzleFuel * rNozzleFuel) / (NUM_OF_TAN_INPUT_NOZZLE_FUEL * A_NOZZLE_FUEL)) = sqrt((%s * %s) / (%s * %s)) = %s мм\n", decimal2Format.format(rInputNozzleFuel * 1000), decimal2Format.format(rNozzleFuel * 1000), NUM_OF_TAN_INPUT_NOZZLE_FUEL, decimal2Format.format(A_NOZZLE_FUEL), decimal2Format.format(rTanInputNozzleFuel * 1000));

        //double reInputNozzleFuel = 4 * mFlowNozzleFuel / ()

        double squareNozzleGG = mFlowNozzleGG / (NU_NOZZLE_GG * ggAstraRes.getRho() * Math.pow(pKs / (pKs + DELTA_P_NOZZLE_GG), 1 / ggAstraRes.getK()) * Math.sqrt((2 * ggAstraRes.getK() / (ggAstraRes.getK() - 1) * ggAstraRes.getT() * ggAstraRes.getK() * (1 - Math.pow(pKs / (pKs + DELTA_P_NOZZLE_GG), (ggAstraRes.getK() - 1) / ggAstraRes.getK())))));
        System.out.printf("squareNozzleGG = %s * 10^-5 м^2\n", decimal2Format.format(squareNozzleGG * 100000));
        double dNozzleGG = Math.sqrt(4 * squareNozzleGG / Math.PI);
        System.out.printf("dNozzleGG = sqrt(4 * squareNozzleGG / PI) = sqrt(4 * %s * 10^-5 / PI) = %s мм\n", decimal3Format.format(squareNozzleGG * Math.pow(10, 5)), decimal2Format.format(dNozzleGG * 1000));
        IevleevResult ievleevResult = getKMByIevlevMethod(mixedHead, dKs, mFlowNozzlePr, mFlowOxInCenter / numOfNozzleOx, mFlowFuelGGCenter / numOfNozzleOx, mFlowNozzleFuel, true);

    }

    private static IevleevResult getKMByIevlevMethod(MixedHead mixedHead, double dKs, double mFlowPr, double mFlowGGOx, double mFlowGGFuel, double mFlowFuel, boolean print) {
        int numOfSteps = 5;
        List<List<Double>> result = new LinkedList<>();

        double rOfArea = mixedHead.getH() * 3;
        double rMax = (dKs - mixedHead.getH()) / 2;
        double rStep = rMax / (numOfSteps - 1);
        double argStep = (Math.PI / 4) / (numOfSteps);
        double arg = 0;
        while (arg - (Math.PI / 4) < 0.0001) {
            double r = 0;
            List<Double> resultByArg = new LinkedList<>();
            while (r <= rMax) {
                Vector<Double> currCenter = new Vector<>(r * Math.cos(arg), r * Math.sin(arg));
                double finalR = r;
                double finalArg = arg;
                List<Vector<Double>> pr = mixedHead.getPrNozzles().stream().filter(v -> checkPointInArea(currCenter, v, rOfArea))
                        .map(v -> toNewCoordSystemInFirstQuarter(v, finalR, finalArg)).toList();
                List<Vector<Double>> fuel = mixedHead.getFuelNozzles().stream().filter(v -> checkPointInArea(currCenter, v, rOfArea))
                        .map(v -> toNewCoordSystemInFirstQuarter(v, finalR, finalArg)).toList();
                List<Vector<Double>> oxAndFuel = mixedHead.getOxNozzles().stream().filter(v -> checkPointInArea(currCenter, v, rOfArea))
                        .map(v -> toNewCoordSystemInFirstQuarter(v, finalR, finalArg)).toList();
                boolean lastStepPerR = r + rStep >= rMax;
                double denominator = mFlowPr * pr.stream().map(v -> getFrenelResult(v, mixedHead.getH(), lastStepPerR)).reduce(Double::sum).orElse(0.0);
                denominator += mFlowFuel * fuel.stream().map(v -> getFrenelResult(v, mixedHead.getH(), lastStepPerR)).reduce(Double::sum).orElse(0.0);
                denominator += mFlowGGFuel * oxAndFuel.stream().map(v -> getFrenelResult(v, mixedHead.getH(), lastStepPerR)).reduce(Double::sum).orElse(0.0);
                double numerator = mFlowGGOx * oxAndFuel.stream().map(v -> getFrenelResult(v, mixedHead.getH(), lastStepPerR)).reduce(Double::sum).orElse(0.0);
                double kM = numerator / denominator;
                resultByArg.add(kM);
                r += rStep;
            }
            result.add(resultByArg);
            arg += argStep;
        }

        if (print) {
            arg = 0;
            for (int i = 0; i < result.size(); i++) {
                double r = 0;
                List<Double> currKmByR = result.get(i);
                for (int j = 0; j < currKmByR.size(); j++) {
                    String formattedX = String.format("%.5f", r).replace(",", ".");
                    String formattedY = String.format("%.5f", currKmByR.get(j)).replace(",", ".");
                    System.out.printf("(%s, %s),", formattedX, formattedY);
                    r += rStep;
                }
                System.out.println("");
                arg += argStep;
            }
        }

        return new IevleevResult(result, argStep, rStep, (Math.PI / 4), rMax);
    }

    private static double getFrenelResult(Vector<Double> point, double h, boolean lastPerR) {
        double hSqrt2 = h * Math.sqrt(2);
        double hDiv2 = h / 2;

        double zx1 = (point.getX() - hDiv2) / hSqrt2;
        double zx2 = (point.getX() + hDiv2) / hSqrt2;
        double zy1 = (point.getY() - hDiv2) / hSqrt2;
        double zy2 = (point.getY() + hDiv2) / hSqrt2;

        return INTEGRATOR.integrate(10000, FRENEL_UNDER_INTEGRATE_FUNC::apply, Math.min(zx1, zx2), Math.max(zx1, zx2))
                * (lastPerR ? (Math.sqrt(Math.PI) / 2) - INTEGRATOR.integrate(10000, FRENEL_UNDER_INTEGRATE_FUNC::apply, 0, Math.min(Math.abs(zy1), Math.abs(zy2))) : INTEGRATOR.integrate(10000, FRENEL_UNDER_INTEGRATE_FUNC::apply, Math.min(zy1, zy2), Math.max(zy1, zy2)))
                * 2 / Math.sqrt(Math.PI);
    }

    private static boolean checkPointInArea(Vector<Double> centerPoint, Vector<Double> point, double r) {
        double x = centerPoint.getX() - point.getX();
        double y = centerPoint.getY() - point.getY();
        return Math.sqrt(x * x + y * y) <= r;
    }

    private static Vector<Double> toNewCoordSystemInFirstQuarter(Vector<Double> point, double r, double arg) {
        double nineteenMinusArg = Math.min((Math.PI / 2) - arg, (Math.PI / 2));
        Vector<Double> newAxiX = new Vector<>(Math.cos(nineteenMinusArg), (-1) * Math.sin(nineteenMinusArg));
        Vector<Double> newAxiY = new Vector<>(Math.cos(arg), Math.sin(arg));
        Vector<Double> res = new Vector<>(point.getX(), point.getY());
        res.setX(res.getX() - r * Math.cos(arg));
        res.setY(res.getY() - r * Math.sin(arg));
        double newX = res.getX() * newAxiX.getX() + res.getY() * newAxiX.getY();
        double newY = res.getX() * newAxiY.getX() + res.getY() * newAxiY.getY();
        res.setX(newX);
        res.setY(newY);
        return res;
    }

    private static MixedHead getMixedHead(
            int nPrSlice, int nPrNozzlePerSlice, int nFuelSlice, double dPr, double dKs, boolean printPoints
    ) {
        List<Vector<Double>> pointsPr = new ArrayList<>();
        double[] dPrSlice = new double[nPrSlice];
        for (int i = 0; i < dPrSlice.length; i++) {
            dPrSlice[i] = dKs - (i + 1) * 2 * (dPr) * COEF_STEP_PR;
            double arg = 0;
            double stepAgr = 2 * Math.PI / nPrNozzlePerSlice;
            while (arg < 2 * Math.PI) {
                pointsPr.add(new Vector<>(dPrSlice[i] * Math.cos(arg) / 2, dPrSlice[i] * Math.sin(arg) / 2));
                arg += stepAgr;
            }
        }
        int numCenterSlice = 2 * nFuelSlice - 1;
        double d = dPrSlice[dPrSlice.length - 1] - 2 * dPr * COEF_STEP_PR;
        double h = d / (2 * (numCenterSlice - 1));
        System.out.printf("H: %s\n", h);
        System.out.printf("d: %s\n", d);
        List<Vector<Double>> pointsFuel = new ArrayList<>();
        List<Vector<Double>> pointsOx = new ArrayList<>();

        double yStep = h * Math.sqrt(3) / 2;
        int yNumSlice = 0;
        double y = 0;
        while (y <= d / 2) {
            double x = yNumSlice % 2 == 0 ? 0 : h / 2;
            int xNumSlice = 0;
            while (Math.sqrt(x * x + y * y) - d / 2 < 0.00001) {
                Vector<Double> point = new Vector<>(x, y);
                List<Vector<Double>> pointsSet;
                if (yNumSlice % 2 == 0) {
                    pointsSet = xNumSlice % 3 == 0 ? pointsFuel : pointsOx;
                } else {
                    pointsSet = (xNumSlice - 1) % 3 == 0 ? pointsFuel : pointsOx;
                }
                pointsSet.add(point);
                if (point.getX() != 0) {
                    pointsSet.add(new Vector<>((-1) * point.getX(), point.getY()));
                }
                if (point.getX() != 0 && point.getY() != 0) {
                    pointsSet.add(new Vector<>((-1) * point.getX(), (-1) * point.getY()));
                }
                if (point.getY() != 0) {
                    pointsSet.add(new Vector<>(point.getX(), (-1) * point.getY()));
                }
                x += h;
                xNumSlice++;
            }
            y += yStep;
            yNumSlice++;
        }
        if (printPoints) {
            List<Vector<Double>> listForPrint = new ArrayList<>(pointsOx);
            listForPrint.addAll(pointsFuel);
            listForPrint.addAll(pointsPr);
            System.out.print("OX:");
            pointsOx.forEach((v) -> {
                String formattedX = String.format("%.5f", v.getX()).replace(",", ".");
                String formattedY = String.format("%.5f", v.getY()).replace(",", ".");
                System.out.printf("(%s, %s),", formattedX, formattedY);
            });
            System.out.print("\nPR:");
            pointsPr.forEach((v) -> {
                String formattedX = String.format("%.5f", v.getX()).replace(",", ".");
                String formattedY = String.format("%.5f", v.getY()).replace(",", ".");
                System.out.printf("(%s, %s),", formattedX, formattedY);
            });
            System.out.print("\nFUEL:");
            pointsFuel.forEach((v) -> {
                String formattedX = String.format("%.5f", v.getX()).replace(",", ".");
                String formattedY = String.format("%.5f", v.getY()).replace(",", ".");
                System.out.printf("(%s, %s),", formattedX, formattedY);
            });
            System.out.println("");
        }
        return new MixedHead(pointsOx, pointsFuel, pointsPr, h);
    }

    private static double getTNABalanceResult(AstraResult astraResultForGG, double mOx, double mFuel, double pKs) {
        System.out.printf("p1O = %s МПа (давление на входе в насос ок) \n", P1_O);
        System.out.printf("p1G = %s МПа (давление на входе в насос гор)\n", P1_G);
        System.out.printf("etaOx = %s (КПД насоса ок)\n", ETA_OX);
        System.out.printf("etaFuel = %s (КПД насоса гор)\n", ETA_FUEL);
        System.out.printf("etaT = %s (КПД турбины)\n", ETA_T);
        System.out.printf("deltaPMagOx = %s МПа (потери на магистрали окислителя)\n", DELTA_P_MAG_OX);
        System.out.printf("deltaPMagFuel = %s МПа (потери на магистрали горючего и форсунки ГГ)\n", DELTA_P_MAG_FUEL);
        System.out.printf("deltaPMagFuel1 = %s МПа (сопротивление магистрали горючего  от выхода из насоса до входа в камеру сгорания (с учетом сопротивления форсунок))\n", DELTA_P_MAG_FUEL1);
        System.out.printf("deltaPT1 = %s МПа (потери на отвод газа турбины к кам сгорания)\n", DELTA_P_T1);
        System.out.printf("deltaPT2 = %s МПа (потери на ответ из ГГ в турбину)\n", DELTA_P_T2);
        double r = astraResultForGG.getR();
        double t = astraResultForGG.getT();
        double k = astraResultForGG.getK();
        double kGG = GG_ALPHA * KM0;

        Function<Double, Double> requiredPower = pGg -> Math.pow(10, 3) * (mOx / (ETA_OX * RHO_OX)) * (pGg + DELTA_P_MAG_OX - P1_O) + (mFuel / (ETA_FUEL * RHO_FUEL)) * (pGg + DELTA_P_MAG_FUEL - P1_G); //кВт
        Function<Double, Double> turbinePower = pGg -> Math.pow(10, -3) * r * t * (1 - Math.pow((pKs + DELTA_P_T1) / (pGg - DELTA_P_T2), (k - 1) / k)) * mOx * ((1 + kGG) / kGG) * ETA_T * k / (k - 1); //кВт

        BrentSolver solver = new BrentSolver(0.01);
        double pGG = solver.solve(100, p -> requiredPower.apply(p) - turbinePower.apply(p), 1, 30);
        System.out.printf("pGG = %s МПа (подбор корней)\n", decimal1Format.format(pGG));

        double p2O = pGG + DELTA_P_MAG_OX;
        double p2G = Math.max(pGG + DELTA_P_MAG_FUEL, pKs + DELTA_P_MAG_FUEL1);
        double p1T = pGG - DELTA_P_T2;
        double p2T = pKs + DELTA_P_T1;
        System.out.printf("p2O = pGG + deltaPMagOx = %s + %s = %s МПа\n", decimal1Format.format(pGG), decimal1Format.format(DELTA_P_MAG_OX), decimal1Format.format(p2O));
        System.out.printf("p2G = max(pGG + deltaPMagFuel; pKs + deltaPMagFuel1) = max(%s + %s; %s + %s) = %s МПа\n",
                decimal1Format.format(pGG),
                decimal1Format.format(DELTA_P_MAG_FUEL),
                decimal1Format.format(pKs),
                decimal1Format.format(DELTA_P_MAG_FUEL1),
                decimal1Format.format(p2G));
        System.out.printf("p1T = pGG - deltaPT2 = %s - %s = %s МПа\n", decimal1Format.format(pGG), decimal1Format.format(DELTA_P_T2), decimal1Format.format(p1T));
        System.out.printf("p2T = pKs + deltaPT1 = %s + %s = %s МПа\n", decimal1Format.format(pKs), decimal1Format.format(DELTA_P_T1), decimal1Format.format(p2T));
        double piT = p1T / p2T;
        System.out.printf("piT = p1T / p2T = %s / %s = %s\n", decimal1Format.format(p1T), decimal1Format.format(p2T), decimal2Format.format(piT));
        return pGG;
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
            if (i == squares.length - 1) {
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

    private static class MixedHead {
        private List<Vector<Double>> oxNozzles;
        private List<Vector<Double>> fuelNozzles;
        private List<Vector<Double>> prNozzles;
        private double h;

        public MixedHead(List<Vector<Double>> oxNozzles, List<Vector<Double>> fuelNozzles, List<Vector<Double>> prNozzles, double h) {
            this.oxNozzles = oxNozzles;
            this.fuelNozzles = fuelNozzles;
            this.prNozzles = prNozzles;
            this.h = h;
        }

        public List<Vector<Double>> getPrNozzles() {
            return prNozzles;
        }

        public void setPrNozzles(List<Vector<Double>> prNozzles) {
            this.prNozzles = prNozzles;
        }

        public List<Vector<Double>> getOxNozzles() {
            return oxNozzles;
        }

        public void setOxNozzles(List<Vector<Double>> oxNozzles) {
            this.oxNozzles = oxNozzles;
        }

        public List<Vector<Double>> getFuelNozzles() {
            return fuelNozzles;
        }

        public void setFuelNozzles(List<Vector<Double>> fuelNozzles) {
            this.fuelNozzles = fuelNozzles;
        }

        public double getH() {
            return h;
        }

        public void setH(double h) {
            this.h = h;
        }
    }

    private static class IevleevResult {
        private List<List<Double>> kms;
        private double argStep;
        private double rStep;
        private double argMax;
        private double rMax;

        public IevleevResult(List<List<Double>> kms, double argStep, double rStep, double argMax, double rMax) {
            this.kms = kms;
            this.argStep = argStep;
            this.rStep = rStep;
            this.argMax = argMax;
            this.rMax = rMax;
        }

        public List<List<Double>> getKms() {
            return kms;
        }

        public void setKms(List<List<Double>> kms) {
            this.kms = kms;
        }

        public double getArgStep() {
            return argStep;
        }

        public void setArgStep(double argStep) {
            this.argStep = argStep;
        }

        public double getrStep() {
            return rStep;
        }

        public void setrStep(double rStep) {
            this.rStep = rStep;
        }

        public double getArgMax() {
            return argMax;
        }

        public void setArgMax(double argMax) {
            this.argMax = argMax;
        }

        public double getrMax() {
            return rMax;
        }

        public void setrMax(double rMax) {
            this.rMax = rMax;
        }
    }

}
