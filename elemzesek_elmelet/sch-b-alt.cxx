/* A Schoenflies-Bieberbach tetelt altalanasithatjuk tetszoleges geometriara
 * ugy, hogy racsok helyett fixpont mentes egybevagosagokat vizsgalunk.
 * 
 * Irunk egy modszert, ami egy adott D-szimbolumhoz eloallit egy fixpont mentes
 * egybevagosagcsoportot leiro "felfujt" D-szimbolumot.
 *
 * Az algoritmus:
 * - Sorra vesszuk a tukrozeseket (lexikografikusan elso az operacio szama,
 * masodik a szimplex sorszama,) es minden tukrozesre megduplazzuk a
 * D-szimbolumot. Az eredeti es az uj D-szimbolumot az adott tukrozeshez tartozo
 * operacio koti csak ossze. (Ilyenkor el tud romlani a nem szomszedos elek
 * korberakasanak feltetele... lasd rajz1) Ez nem eleg.
 * - Javitott elso lepes: Megnezzuk az ekvivalens tukrozeseket. Mit is jelent
 * ez?
 *   * a sulyzo alaku rendszereket ki kell javitani vele (rajz2)
 *
 * - Uj verzio: Az osszes elet nezzuk:
 *   * Sorra vesszuk a nem iranyitottan (tukrozes segitsegevel) korberakott
 *   eleket. Minden ilyet megszuntetunk ugy, hogy a korberakashoz tartozo
 *   parameter ketszereseszer lemasoljuk a teljes D-szimbolumot, es a tukrozott
 *   szimplexeket es masolataikat kapcsoljuk ossze a tukrozes operatoraval
 *   (rajz3) Ez igy nem jo!!! El tudja rontani a nem szomszedos operator-parok
 *   rendjet. Mi a megoldas? 
 *     o Kivalasztjuk az egyik szimplexre hato egyik tukrozes operatort (amit
 *     meg akarunk szuntetni)
 *     o Megnezzuk, hogy ez az operator melyik elek korberakasaban jatszik
 *     szerepet.
 *     o Mikor romlik el a helyzet?
 *       x sulyzo alak (nem szomszedos operaciokkal; vagy szomszedos, de 6-nal
 *       kisebb M ertekkel)
 *       x Egy hurkot szuntetunk meg, ezert minden hozza kapcsolodo
 *       matrix-ertekhez vagy egy (valamilyen hosszu) sulyzo, vagy egy
 *       dupla-hurok tartozik
 *     o A megoldas:
 *       x Az M/R ertekeket figyeljuk, amennyiben ez 1 valamelyik operacio-par
 *       eseten, akkor a teljes kort bezarjuk...
 *     o Ujabb problema: Mi van, ha nem paros a parameter (mar csak a szomszedos
 *     operaciok johetnek szoba...) Ekkor ugyanis nem erunk vissza az eredeti
 *     szimplex-osztalyba.
 *       x Durvabb hozaallas: Kivalasztunk egy tukrozest es minden elt
 *       korberakunk, amiben szerepel. Hat nem...
 *       x A helyzet meg rosszabb, ha a 3-3-3 parameteru, 1 elemu D-szimbolumot
 *       nezzuk (Minden parameter paratlan, ha az egyiket befoltozom, a masik
 *       elkezd ereszteni...) Ez egy gombi csoportot hataroz meg, ami pedig
 *       veges, tehat eloszor fel kell ismerni a geometriat, kulonben nem tudjuk
 *       megmondani, hogy az ismetlodes miert van...
 *     o Kene talalni kicsi lepeseket, mert az egeszet nem fogjuk tudni
 *     felismerni egy lepesben. Kis lepesekhez valamilyen allandosagra van
 *     szukseg. Ezt nem talalom... Lehet, hogy a fundamentalis tartomanyra lesz
 *     szuksegem, hiszen tukrozeseket, forgatasokat akarok eltuntetni.
 *       x Kicsi lepesek, amik nem rontjak el a parametereket:
 *       Mindig egy nem szomszedos operacio-parhoz tartozo sulyzot (vagy dupla
 *       hurkot) bontunk szet. Ezzel megtobbszorozzuk (vagy megparos-szorozzuk) a
 *       D-szimbolumot. Ezeket a lepeseket addig folytatjuk, amig a megfelelo
 *       parametert el nem erjuk, ekkor visszaterunk a kiindulasi ponthoz. (Ez
 *       akkor jo, ha 2 oldala van a szimbolumnak (lasd abra5.)
 *     o Egy nagyobb D-szimbolumban egyszeruen eszre lehet venni egy kisebbet
 *     (lasd maximalis-e kerdes a diplomamunkaban.) Ezt nem lehet megforditani?
 *   * Sorra vesszuk az iranyitottan (tukrozes nelkul) korberakott eleket.
 *   Minden ilyet megszuntetunk ugy, hogy eloszor szetvagjuk a kort egy
 *   tetszoleges helyen, majd a korberakashoz tartozo parameterszer lemasoljuk a
 *   teljes szetvagot D-szimbolumot, es a vagashoz betoldjuk (rajz4)
 */
