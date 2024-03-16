package ru.bmstu.krdu.project.dto;

public class Vector<Y> {
    private double x;
    private Y y;

    public Vector(double x, Y y) {
        this.x = x;
        this.y = y;
    }

    public double getX() {
        return x;
    }

    public void setX(double x) {
        this.x = x;
    }

    public Y getY() {
        return y;
    }

    public void setY(Y y) {
        this.y = y;
    }

}
