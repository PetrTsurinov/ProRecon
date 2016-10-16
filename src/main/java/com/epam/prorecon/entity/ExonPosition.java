package com.epam.prorecon.entity;

public class ExonPosition {
    private String chromosome;
    private ExonType exonType;
    private int startIndex;
    private int endIndex;

    public ExonPosition(String chromosome, ExonType exonType, int startIndex, int endIndex) {
        this.chromosome = chromosome;
        this.exonType = exonType;
        this.startIndex = startIndex;
        this.endIndex = endIndex;
    }

    public ExonPosition() {
    }

    public ExonPosition(ExonPosition exonPosition) {
        this.chromosome = exonPosition.chromosome;
        this.exonType = exonPosition.exonType;
        this.startIndex = exonPosition.startIndex;
        this.endIndex = exonPosition.endIndex;
    }

    public void setExonPosition(ExonPosition exonPosition) {
        this.chromosome = exonPosition.chromosome;
        this.exonType = exonPosition.exonType;
        this.startIndex = exonPosition.startIndex;
        this.endIndex = exonPosition.endIndex;
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + (chromosome != null ? chromosome.hashCode() : 0);
        result = prime * result + (exonType != null ? exonType.hashCode() : 0);
        result = prime * result + new Integer(startIndex).hashCode();
        result = prime * result + new Integer(endIndex).hashCode();
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        ExonPosition other = (ExonPosition) obj;

        return ((chromosome == null && other.chromosome == null) || (chromosome != null && chromosome.equals(other.chromosome)))
                && ((exonType == null && other.exonType == null) || (exonType != null && exonType.equals(other.exonType)))
                && startIndex == other.startIndex && endIndex == other.endIndex;
    }

    public ExonType getExonType() {
        return exonType;
    }

    public int getStartIndex() {
        return startIndex;
    }

    public void setStartIndex(int startIndex) {
        this.startIndex = startIndex;
    }

    public int getEndIndex() {
        return endIndex;
    }

    public void setEndIndex(int endIndex) {
        this.endIndex = endIndex;
    }
}
