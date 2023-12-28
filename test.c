#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>
#include <stdarg.h>
#include "LLL_Reduction.h"
#include "Enumeration.h"
#include <stdbool.h>


double SumArray(int dim, double arr[dim]) {
    int i;
    double sum1;
    for (i=0; i<dim; i++) {
        sum1 += arr[i];
    }
    return sum1;
}

double Determinant(int dim, double (*A)[dim]) {
    int i, j, k;
    double sum1;
    bool skip[2] = {false, false};
    double B[dim-1][dim-1];
    if (dim == 1) {
        return A[0][0];
    }
    for (i=0; i<dim; i++) {
        skip[0] = skip[1] = false;
        for (j=0; j<dim; j++) {
            if (j==i) {
                skip[0] = true;
            }
            for (k=0; k<dim; k++) {
                if (k==i) {
                    skip[1] = true;
                }
                B[j][k] = A[j+skip[0]][k+skip[1]];
            }
        }
        sum1 += pow(-1, i) * A[0][i] * Determinant(dim-1, B);
    }
}

double LimitCalc(int dim, double (*A)[dim]) {
    double gamma = tgamma(dim/2 + 1);
    double det = Determinant(dim, A);
    return 1.05*(pow(gamma, 1/dim)/sqrt(M_PI))*pow(det, 1/dim);
}

void runTests(int dim, ...)
{
    //typedef struct {
    //    double *array;
    //    size_t size;
    //} Array;
    
    int i, j;
    double A[dim][dim];
    double *arr;
    double limit;
    
    va_list args;
    va_start (args, dim);
    for (i = 0; i < dim; i++) {
        arr = va_arg(args, double*);
        for (j=0; j < dim; j++) {
            A[i][j] = arr[j];
        }
    }
    printf("A\n");
    for (i = 0; i < dim; i++) {
        for (j=0; j < dim; j++) {
          printf("%.4f\t", A[i][j]);
        }
        printf("\n");
    }
    
    LLL(0.75, 3, A);
    printf("Orthonormalized Vectors (A):\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%.4f\t", A[i][j]);
        }
        printf("\n");
    }
    double shortest_vector = ShortestVector(3, A);

    bool unit_test = true;
    for (i=0; i<dim; i++) {
        if (SumArray(dim, A[i]) != 1.0) {
            unit_test = false;
            break;
        }
    }

    if (unit_test) {
        printf("For Dimension: %d Expected: %.4f Got: %.4f\n", dim, 1.0, shortest_vector);
        if(shortest_vector != 1.0) {
            printf("Expected %.4f, got %.4f\n", 1.0, shortest_vector);
        }
        assert(shortest_vector == 1.0);
    }

    else {
        limit = LimitCalc(dim, A);
        printf("For Dimension: %d Limit: %.4f Got: %.4f\n", dim, limit, shortest_vector);
        if(shortest_vector > limit) {
            printf("Limit %.4f, got %.4f\n", limit, shortest_vector);
        }
        assert(shortest_vector <= limit);
    }
}

int main() {

    
    runTests(40, {371644562438531748585630085667746983470408244373580742281269821738816154869236300428137729771977309321210236491639933332 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {970128588937842140661210199069589679615703684562655262446230094731888177573698872897660975186881893815551083050469479861 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {1617885213873858163065250499756376305837976609435400008057297118650398723827143214937143752380433141030895046085134663139 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {861650440856876922136315403300972806928844583476323346957866774048720276102084753290592827764541119558381900492846561687 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {1238401879408921089844474386903737818276239209740348854610933282304058072380604156331847453869318599590680629601714021166 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {302328193892196204243023116106370462383132199983296195152211088133814423239197999286069616594274614510736145263925517783 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {443681408734567838343607093614379754005413471877157728597635239760454843085205167872415159913459109353288267514781828878 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {766924867740885236236436059448427698090379772276332984401846634621581522685296543079817924330812470098383411173196131706 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {1297358997739527278188204320802855592592899794686777403497982861534586993960747704855480346961487202646587358008424555998 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {1923256447555821027760823712014018288715572436486045654757604016247481734877091777864720799782922934671078042956616750575 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {138791684593421316484912136928975449358519306209923545534576080943933941788220318189376736295439954669622762228268690319 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {1752053685015713567891789191640042337741805857685242642164798009587123835789029180360779591122248108173121788944806990018 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {2046303652643607551083773382245919660376167541225343779851326688493989336785866602067814167781776230842039502024961224990 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {161701314003903527294539093851830420185672588693938679547737547227941952578846685981290488906797666486385636306427956922 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {2059714413420462412604210521069862437228594984719691499291824764864200097906908418575923366701593028750494787160369423664 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {236619910841436753359639417128223769093980125208659184950079289394351446700013212860226717401403083201476790108710198995 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {521807735461695034807439230274452751978618453371762238147892673540779228084610103951896314431773738642862361617950116793 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {1898466106725110481310768025369914138166711138228177956551055744552328440573249624614143156799333166991080740657486025033 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {2043257177335956951070219397858777123559357103945545237417070888207682821645422085875324039671029805011652245597869305538 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {904225163671074382768019753539274681795779625382511961264465653291127670514453520050551644296901230356912243654225654401 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {32771872907828211382289261014117057923327418144915569087238809343521061214616361304729214156019094728432041297468304493 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {185789413730166201171458176468755969088909804647627433894406398430754979813498091666025324336151794710682026172503101451 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {832962803741579150035009987801883256172874768041791705764900197393212435979706525525238578679210561346645899576161141793 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {1042621519072433920283083562428469593602182919945835211322354565330887090220897472063018488024892773750305996602609483380 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {171388995082154233233012548090942599143715724202696890942088990082340274801161970138190639906868753248272491463187401695 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0}, {2100493671200104536533474921602153259773348308608862002579228988953071977809598351863259953863461834410595513598037485038 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0}, {849028874593296703934398545250196283716795306401228544442045159484861904205470885973430411009159267090428932112495507852 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0}, {1046701824957941065045627030825750373570344803958334982369286898950041088016988929016615138537061510672807965422152953094 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0}, {542875223055107687983299240475304624160994620138110511564200316784444891282138558584574906974554732933893987639864801905 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0}, {1584382682441666539891259939997296576922760554674178579888435692973226441994929350838648607785926118636355860195320578386 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0}, {316503555252451563889967315570982477067169012390726271667771540315451791264064871351402596704077496942199692552980634552 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0}, {1837738370141979294787389742951594377674155679679663044021277108096701623188534975556948353238304199946688862260850711361 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0}, {739578827111769709363025786954762546314747907261948920241963639022818938531742051897853101539916787244178755946187647555 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0}, {1619290331242112257658768582175275050536918255267445777163956250682374601071625719555758108294380161274992932056486622464 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0}, {852244520379043286720418573598354118899465075623425806989984806937105488837636350325384144744195783077780523735920132437 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0}, {491483555886899765880626820064070785891162464132735054675483057099729053163156386346862511333336434543264870015359686516 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0}, {842211543902380452313476978966614900223458281660135345978369910918684855621079286840261544398613444576199814061553688778 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0}, {1874411128584029818874533765431554481569149176662206972692659544113727361200950942823993552649265946057735986943133343384 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0}, {1807564784681621941238943795385610371403525911952667414563776455384585295644151678837176854600997116364524617338835884611 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1});
    
    //runTests(7, 216);
    //runTests(20, 114624);
    //runTests(30, 14098308);
    return 0;
}
