// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace std;
using namespace Eigen;


// [[Rcpp::export(name="BODTCpp")]]
NumericVector GBDT(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::MatrixXd Xtest, NumericMatrix Y,
                   int numClass = 1, int maxTerms = 100, int ntrees = 100, int mtry=1, int MinLeaf = 5,
                   bool replacement=true, double ratOOB=0.37)
                   //Rcpp::Function nnet, Rcpp::Function rpart, Rcpp::Function control, Rcpp::Function predict,
{

  int i, j, ii, jj, k, t, nterm, nterm_fails = 0, n = y.size(), nn=n, n1 = Xtest.rows(), p = X.cols();
  double mse, aic, AIC = 1e+10;
  if(!replacement)nn=ceil((1-ratOOB)*n);

  Rcpp::IntegerVector TDindx, J, NTD, sn = Rcpp::seq(0, n - 1), sp = Rcpp::seq(0, p - 1);
  Eigen::MatrixXd XBtest(n1, mtry), XBsub(nn, mtry), XB(n, mtry), B(mtry, 1);

  Rcpp::Environment pkg1=Environment::namespace_env("nnet"), pkg2=Environment::namespace_env("rpart");// pkg3=Environment::namespace_env("stats");
  Rcpp::Function nnet=pkg1["nnet.default"], control=pkg2["rpart.control"], rpart=pkg2["rpart"],predict=pkg2["predict.rpart"];// lm=pkg3["lm"];
  Rcpp::Formula formula = Rcpp::Formula("y~.");
  Rcpp::List fit;

  if (numClass > 2)
  {
    nterm = 0;
    NumericMatrix YRes(n, numClass),YBRes(nn, numClass), RES(n, numClass),COEF(maxTerms + 1, numClass);
    Rcpp::List fitted(numClass), pred(numClass), XXBsub,XXB, XXBtest;                    // create a list with numMatrices elements
    NumericMatrix fittedC(n, maxTerms), predC(n1, maxTerms), COEF0; // create a matrix
    NumericVector coef(maxTerms + 1);

    Rcpp::Environment pkg3 = Environment::namespace_env("stats");
    Rcpp::Function lm = pkg3["lm"];

    // fill the list with copies of the matrix
    for (i = 0; i < numClass; i++)
    {
      fitted[i] = clone(fittedC); // use clone() to make a copy of the matrix
      pred[i] = clone(predC);
    }

    YRes = clone(Y);//;//Y; //

    for (t = 0; t < maxTerms; t++)
    {

      if (ntrees == 1)
      {
        TDindx = sn; // Rcpp::clone(sn);
        J = sp;
      }
      else
      {
        J = Rcpp::sample(sp, mtry, FALSE);
        TDindx = Rcpp::sample(sn, nn, replacement);
      }
      NTD = Rcpp::setdiff(sn, TDindx);

      jj = 0;
      std::for_each(std::begin(J), std::end(J), [&](int j){
        XB.col(jj)=X.col(j);
        XBtest.col(jj)=Xtest.col(j);
        jj++; });

      ii = 0;
      std::for_each(std::begin(TDindx), std::end(TDindx), [&](int i){
        XBsub.row(ii)=XB.row(i);
        YBRes.row(ii)=YRes.row(i);
        ii++; });

      fit = nnet(XBsub, Rcpp::Named("y", YBRes), Rcpp::Named("size", 1), Rcpp::Named("trace", false), Rcpp::Named("linout", true)); // linout = TRUE,MaxNWts= p+5
      B.col(0) = (as<Eigen::VectorXd>(fit["wts"])).segment(1, mtry);                                                            //.head(mtry);
      XXBsub = DataFrame::create(XBsub, _["XB"] = XBsub * B);
      XXB = DataFrame::create(XB, _["XB"] = XB * B);
      XXBtest = DataFrame::create(XBtest, _["XB"] = XBtest * B);

      k=0;
      for (k = 0; k < numClass; k++)
      {
        fit = rpart(formula, Rcpp::Named("data") = DataFrame::create(XXBsub, _["y"] = YBRes.column(k)),
                    Rcpp::Named("control") = control(Rcpp::Named("minbucket", MinLeaf)));

        // NumericMatrix fittedk=fitted[k];
        fittedC = wrap(fitted[k]);//wrap();
        fittedC.column(nterm) = as<NumericVector>(predict(fit, Rcpp::Named("newdata") =XXB));
        // fitted[k]=fittedC;
        predC = wrap(pred[k]);
        predC.column(nterm) = as<NumericVector>(predict(fit, Rcpp::Named("newdata") = XXBtest));

        fit = lm(formula, Rcpp::Named("data") = DataFrame::create(fittedC(_, Range(0, nterm)), _["y"] = Y.column(k)));
        RES.column(k) = as<NumericVector>(fit["residuals"]);

        // fit["coefficients"][Rcpp::is_na(fit["coefficients"])] <- 0.0;
        // COEF.column(k)[Range(0,nterm)] = fit["coefficients"];
        //coef = COEF.column(k);
        //COEF(Range(0, nterm+1),k) = as<NumericVector>(fit["coefficients"]);
        coef[Range(0, nterm+1)] = as<NumericVector>(fit["coefficients"]);
        coef[Rcpp::is_na(coef)] = 0.0;
        COEF.column(k)=coef;
      }

      ii = 0;
      aic = sum(RES * RES);

      std::for_each(std::begin(NTD), std::end(NTD), [&](int i)
      { aic += sum(RES.row(i)*RES.row(i)); });

      aic = log(aic / ((n + NTD.size()) * numClass)) + log(p) * nterm * log(n) / n;
      if (aic < AIC)
      {
        YRes = RES;
        AIC = aic;
        COEF0 = COEF(Range(0, nterm+1), _);
        nterm++;
        nterm_fails = 0;
      }
      else
      {
        nterm_fails++;
        if (nterm_fails > 5)
        {
          break;
        }
      }
    }
    nterm--;

    Eigen::VectorXd predictions(n1 * numClass + 1);
    Eigen::MatrixXd PRED(n1, nterm + 2), FITTED(n, nterm + 2), fittes(n, numClass);
    PRED.col(0).setOnes(), FITTED.col(0).setOnes();
    Eigen::Map<Eigen::MatrixXd> BETA(as<Eigen::Map<Eigen::MatrixXd> >(COEF0));
    //Eigen::MatrixXd BETA = as<Eigen::MatrixXd>(COEF0);
    k=0;
    for (k = 0; k < numClass; k++)
    {
      // predC = as<Eigen::Map<Eigen::MatrixXd> >(wrap(pred[k]));
      //predC = wrap(pred[k]);
      //predictions[Range(k * n1, (k + 1) * n1 - 1)] = predC(_, Range(0, nterm))%*%COEF0.column(k);

      PRED.rightCols(nterm+1) = (as<Eigen::Map<Eigen::MatrixXd> >(pred[k])).leftCols(nterm+1);
      //predictions[seq(k * n1, (k + 1) * n1 - 1)] = PRED * BETA.col(k);
      predictions.segment(k * n1, n1) = PRED * BETA.col(k);

      FITTED.rightCols(nterm+1) = (as<Eigen::Map<Eigen::MatrixXd> >(fitted[k])).leftCols(nterm+1);
      fittes.col(k) = FITTED * BETA.col(k);
      //fittedC = wrap(fitted[k]);
      //FITTED.rightCols(nterm) = as<Eigen::MatrixXd>(fittedC(_, Range(0, nterm - 1)));
      //fittes.cols(k) = wrap(FITTED * BETA.col(k));
    }

    fittes=fittes.array()-(as<Eigen::Map<Eigen::MatrixXd> >(Y)).array();
    mse = fittes.array().square().mean();
    //  mse = Rcpp::mean((wrap(fittes) - Y) * (wrap(fittes) - Y));
    predictions(n1 * numClass)=mse;
    //
    return wrap(predictions);


    //Rcpp::List output = Rcpp::List::create(Rcpp::Named("coef")= Y,
    //Rcpp::Named("fitted")=fittes, _["y"] =mse,Rcpp::Named("Y")= predictions); //VectorXi::LinSpaced(n1,k * n1, (k + 1) * n1 - 1)

    //return output;//wrap(B);
  }
  else
  {
    nterm = 1;
    Rcpp::NumericVector res, yres, predictions;
    Eigen::VectorXd betahat, coef;
    Eigen::MatrixXd fitted(n, maxTerms + 1), pred(n1, maxTerms + 1);
    fitted.col(0).setOnes();
    pred.col(0).setOnes();

    if (numClass == 2)
    {
      y = y.array() - 1;
    }
    yres = Rcpp::wrap(y);

    t = 0;
    for (t = 0; t < maxTerms; t++)
    {

      if (ntrees == 1)
      {
        TDindx = sn; // Rcpp::clone(sn);
        J = sp;
      }
      else
      {
        J = Rcpp::sample(sp, mtry, FALSE);
        TDindx = Rcpp::sample(sn, nn, replacement);
      }
      NTD = Rcpp::setdiff(sn, TDindx);


      jj = 0;
      std::for_each(std::begin(J), std::end(J), [&](int j){
        XB.col(jj)=X.col(j);
        XBtest.col(jj)=Xtest.col(j);
        jj++; });

      ii = 0;
      std::for_each(std::begin(TDindx), std::end(TDindx), [&](int i){
        XBsub.row(ii)=XB.row(i);
        ii++; });

      //XB = X(as<vector<int> >(TDindx), as<std::vector<int> >(J));
      // XBtest= XBtest(sn1, J);
      //XB = X(as<std::vector<int> >(J),{4,2,5,5,3});

      //Rcpp::print(X);
      //Rcpp::print(Xtest);
      //Rcpp::Rcout << "x = " << X << std::endl;


      fit = nnet(XBsub, Rcpp::Named("y", yres[TDindx]), Rcpp::Named("size", 1), Rcpp::Named("trace", false), Rcpp::Named("linout", true)); // linout = TRUE,MaxNWts= p+5
      B.col(0) = (as<Eigen::VectorXd>(fit["wts"])).segment(1, mtry);                                                                    //.head(mtry);

      fit = rpart(formula, Rcpp::Named("data") = DataFrame::create(XBsub, _["XB"] = XBsub * B, _["y"] = yres[TDindx]),
                  Rcpp::Named("control") = control(Rcpp::Named("minbucket", MinLeaf)));
      fitted.col(nterm) = as<Eigen::VectorXd>(predict(fit, Rcpp::Named("newdata") = DataFrame::create(XB, _["XB"] = XB * B)));
      pred.col(nterm) = as<Eigen::VectorXd>(predict(fit, Rcpp::Named("newdata") = DataFrame::create(XBtest, _["XB"] = XBtest * B)));

      // fit = lm(formula, Rcpp::Named("data")= DataFrame::create(fitted.leftCols(nterm),_["y"]= y));
      // LM = fastLm(y~.,data.frame(y=yy, fitted))#,silent = TRUE)
      // LM <- ols.fit(x = cbind(1,fitted), y = yy)
      // fit = QRLSCpp(fitted.leftCols(nterm),y);
      // fit= QRLSCpp(fitted.leftCols(nterm+1),y);
      // betahat=fit["coefficients"];
      // res=fit["residuals"];

      HouseholderQR<Eigen::MatrixXd> QR(fitted.leftCols(nterm + 1)); // const
      betahat = QR.solve(y);
      res = y - fitted.leftCols(nterm + 1) * betahat;

      // res=fit["residuals"];
      aic = log((sum(res * res) + sum(res[NTD] * res[NTD])) / (n + NTD.size())) + log(p) * nterm * log(n) / n;
      if (aic < AIC)
      {
        yres = res;
        coef = betahat;
        AIC = aic;
        nterm++;
        nterm_fails = 0;
      }
      else
      {
        nterm_fails++;
        if (nterm_fails > 5)
        {
          break;
        }
      }
    }

    nterm--;
    predictions = pred.leftCols(nterm + 1) * coef;
    mse = (fitted.leftCols(nterm + 1) * coef - y).array().square().mean();
    predictions.push_back(mse);

    return predictions;
  }

}
