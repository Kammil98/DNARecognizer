package dnarecognizer;

public class LevenshteinDistance {
    public static int Compute(String s, String t) {
        if (s.isEmpty()) {
            if (t.isEmpty())
                return 0;
            return t.length();
        }

        if (t.isEmpty()) {
            return s.length();
        }

        int n = s.length();
        int m = t.length();
        int[][] d = new int[n + 1][m + 1];

        // initialize the top and right of the table to 0, 1, 2, ...
        for (int i = 0; i <= n; d[i][0] = i++) ;
        for (int j = 1; j <= m; d[0][j] = j++) ;

        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                int cost = (t.charAt(j - 1) == s.charAt(i - 1)) ? 0 : 1;
                int min1 = d[i - 1][j] + 1;
                int min2 = d[i][j - 1] + 1;
                int min3 = d[i - 1][j - 1] + cost;
                d[i][j] = Math.min(Math.min(min1, min2), min3);
            }
        }
        return d[n][m];
    }
}
