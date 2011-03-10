/* Adott D-szimbolumhoz tartozo fundamentalis tartomany eloallitasa. 
 * 
 * Allitolag megteheto konvex modon. Kerdes, hogy a ragasztas soran vegig konvex
 * marad-e a rendszer, vagy csak a legvegen. (Diplomamunka-beli D.109 eseten pl.
 * nem tud a rendszer vegig konvex maradni.)
 *
 * Mikor tud elromlani a konvexitas (ezeket a helyzeteket ki kell kuszobolni): 
 * - Pontosan akkor lesz az alakzat konkav, ha van konkav lapszoge
 * - Ha adott egy el, amit szimplexekkel korbe tudunk rakni ismetles nelkul, a
 * szimplexek "felenel tobbet", de nem az osszeset ragasztjuk korbe az el korul.
 * Ilyen eleket az iranyitott (tukrozes mentes) 1 erteku paramaterek jelentenek.
 * 
 * Tovabbi gondolatok:
 * - Tukorkepeket nem ragasztunk (ugyanazt a szimplexet 1-szer hasznaljuk fel)
 * - Ha egy el korulrakasaban van tukrozes, nem okoz gondot az el.
 * - Ha nincs tukrozes, de legalabb 2-szer korbemegyunk (parameter min. 2),
 * szinten nem gond.
 * - Ha (akar tukrozessel) 1 elet korbe lehet rakni a szimplexekkel, akkor
 * minden ok (ez konvex lesz).
 * - Azokat az eseteket kell csak megvizsgalni, ahol tobb elet is korbe lehet
 * rakni szimplexekkel ismetles nelkul.
 *
 * Az utolso tovabbi gondolat folytatasa es egy leendo algoritmus:
 * - Vesszuk az iranyitott 1 erteku parametereket (el-korok), az altaluk
 * erintett szimplexek szama szerint csokkeno sorrendben
 * - Az elso ilyenhez tartozo elet korbe ragasztjuk.
 * - Megnezzuk, hogy a maradek el-korok kozul hany olyan van, aminek csak 1
 * kozos szimplexe van az eddig kivalasztott korok uniojaval, vesszuk a
 * legnagyobb ilyet es kezdjuk elolrol az algoritmust (lasd rajz1.)
 * - A maradek el-koroknek legalabb ket szimplexe van az eddig osszeragasztott
 * rendszerben. Vesszuk azt, aminek a legtobb szimplexe van meg szabadon, es
 * ketoldalrol elkezdjuk osszerakni a kort (lasd rajz2.) Fontos, hogy a maradek
 * kort nem tudjuk teljesen osszeragasztani mert csak ekvivalens eleket jar
 * korbe, de nem ugyanazt.
 */
