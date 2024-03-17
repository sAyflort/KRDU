package ru.bmstu.krdu.project.fileUtil;

import ru.bmstu.krdu.project.dto.AstraResult;

import java.text.DecimalFormat;
import java.util.*;

public class AstraUtil {
    public static Map<String, Double> tableOfOptimum(List<AstraResult> astraResults, double foth, double[] alphas, double km0) {
        Map<String, Double> result = new HashMap<>();
        List<AstraResult> ksList = new LinkedList<>();
        List<AstraResult> crList = new LinkedList<>();
        List<AstraResult> outList = new LinkedList<>();
        astraResults.forEach(ar -> {
            if (ar.getFoth() == 0) {
                ksList.add(ar);
            } else if (Math.abs(ar.getFoth() - foth) < 1.5) {
                outList.add(ar);
            } else {
                crList.add(ar);
            }
        });
        if (outList.size() == 1) {
            return null;
        }
        int idx = 0;
        for (int i = 1; i < outList.size(); i++) {
            if (outList.get(i).getI() > outList.get(i - 1).getI()) {
                idx = i;
            }
        }
        DecimalFormat decimalAlphaFormat = new DecimalFormat("0.00");
        DecimalFormat decimalIFormat = new DecimalFormat("0.00");
        System.out.println("OPTIMAL ALPHA");
        System.out.println("KM - ALPHA - I:");
        for (int i = 0; i < outList.size(); i++) {
            System.out.printf("%-6s %-6s %-8s\n", decimalAlphaFormat.format(km0 * alphas[i]), decimalAlphaFormat.format(alphas[i]), decimalIFormat.format(outList.get(i).getI() * 9.81));
        }
        System.out.println("Optimal alpha: " + alphas[idx]);
        System.out.println("-------------------------------------");
        result.put("i", outList.get(idx).getI());
        result.put("fFlowCr", crList.get(idx).getFflow());
        result.put("fFlowOut", outList.get(idx).getFflow());
        result.put("k", ksList.get(idx).getK());
        result.put("pa", outList.get(idx).getP());
        result.put("p", ksList.get(idx).getP());
        result.put("alpha", alphas[idx]);
        return result;
    }

    public static Map<String, Double> tableOfPristenok(List<AstraResult> astraResults, double foth, double[] alphas, double km0) {
        Map<String, Double> result = new HashMap<>();

        List<AstraResult> KSres = new LinkedList<>();
        List<AstraResult> OUTres = new LinkedList<>();

        astraResults.forEach(a -> {
            if (a.getFoth() == 0) {
                KSres.add(a);
            } else if (Math.abs(a.getFoth() - foth) < 1.5) {
                OUTres.add(a);
            }
        });

        int idx = 0;
        for (int i = 1; i < KSres.size(); i++) {
            if (KSres.get(i).getT() > KSres.get(i - 1).getI() && KSres.get(i).getT() <= 2000) {
                idx = i;
            }
        }
        DecimalFormat decimalAlphaFormat = new DecimalFormat("0.000");
        DecimalFormat decimalTFormat = new DecimalFormat("0");
        DecimalFormat decimalIFormat = new DecimalFormat("0.00");
        System.out.println("ALPHA ON WALL");
        System.out.println("KM - ALPHA - T - I:");
        for (int i = 0; i < KSres.size(); i++) {
            System.out.printf(
                    "%-6s %-6s %-6s %-6s\n",
                    decimalAlphaFormat.format(km0 * alphas[i]),
                    decimalAlphaFormat.format(alphas[i]),
                    decimalTFormat.format(KSres.get(i).getT()),
                    decimalIFormat.format(OUTres.get(i).getI() * 9.81)
            );
        }
        System.out.println("Alpha on wall: " + alphas[idx]);
        result.put("t", KSres.get(idx).getT());
        result.put("i", OUTres.get(idx).getI());
        result.put("alpha", alphas[idx]);
        System.out.println("-------------------------------------");
        return result;
    }
}