package ru.bmstu.krdu.project.dto;

public class AstraResult {
    private double p;
    private double k;
    private double t;
    private double r;
    private double i;
    private double foth;

    private double fflow;

    public AstraResult() {
    }

    public double getP() {
        return p;
    }

    public void setP(double p) {
        this.p = p;
    }

    public double getK() {
        return k;
    }

    public void setK(double k) {
        this.k = k;
    }

    public double getT() {
        return t;
    }

    public void setT(double t) {
        this.t = t;
    }

    public double getR() {
        return r;
    }

    public void setR(double r) {
        this.r = r;
    }

    public double getI() {
        return i;
    }

    public void setI(double i) {
        this.i = i;
    }

    public double getFoth() {
        return foth;
    }

    public void setFoth(double foth) {
        this.foth = foth;
    }

    public double getFflow() {
        return fflow;
    }

    public void setFflow(double fflow) {
        this.fflow = fflow;
    }

    public void setField(String name, double value) {
        switch (name) {
            case "p": setP(value); break;
            case "k": setK(value); break;
            case "t": setT(value); break;
            case "r": setR(value); break;
            case "i": setI(value); break;
            case "foth": setFoth(value); break;
            case "fflow": setFflow(value); break;
        }
    }
}
